////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.C,v 1.3 2003-05-23 21:51:04 fgygi Exp $


#include "Species.h"
#include "SpeciesReader.h"
#include "StructuredDocumentHandler.h"
#include "SpeciesHandler.h"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
using namespace std;

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
using namespace xercesc;

////////////////////////////////////////////////////////////////////////////////
SpeciesReader::SpeciesReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::readSpecies (Species& sp, const string uri)
{
  const char* encodingName = "UTF-8";
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  
  int ierr = 0;
  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    cout << "  <!-- Species::readSpecies: Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << " -->" << endl;
    throw;
  }

  parser = XMLReaderFactory::createXMLReader();
  if (valScheme == SAX2XMLReader::Val_Auto)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
      parser->setFeature(XMLUni::fgXercesDynamic, true);
  }

  if (valScheme == SAX2XMLReader::Val_Never)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
  }

  if (valScheme == SAX2XMLReader::Val_Always)
  {
      parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
      parser->setFeature(XMLUni::fgXercesDynamic, false);
  }

  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
  parser->setFeature(XMLUni::fgXercesSchema, doSchema);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
  parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

  int errorCount = 0;
 
  int nlink = 0;
  string current_uri = uri;
 
  while ( sp.description() == "undefined" && nlink++ < 5 )
  {
    try
    {
      SpeciesHandler* sp_handler = new SpeciesHandler(sp);
      StructuredDocumentHandler handler(sp_handler);
      parser->setContentHandler(&handler);
      parser->setErrorHandler(&handler);
      parser->parse(current_uri.c_str());
      errorCount = parser->getErrorCount();
      delete sp_handler;
    }

    catch (const XMLException& toCatch)
    {
        cout << "\nAn error occurred\n  Error: "
             << StrX(toCatch.getMessage())
             << "\n" << endl;
        XMLPlatformUtils::Terminate();
        delete parser;
        throw;
    }

    catch (const SAXParseException& e)
    {
        cout << "\na SAXParseException occurred in file "
             << StrX(e.getSystemId())
             << ", line " << e.getLineNumber()
             << ", char " << e.getColumnNumber()
             << "\n  Message: " << StrX(e.getMessage()) << endl;
        XMLPlatformUtils::Terminate();
        delete parser;
        throw;
    }
 
    if ( sp.description() == "undefined" )
      current_uri = sp.uri();
  }
  sp.uri_ = current_uri;

  delete parser;
  XMLPlatformUtils::Terminate();
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::bcastSpecies(Species& sp)
{
  //cout << ctxt_.mype() << ": starting bcastSpecies" << endl;
  if ( ctxt_.onpe0() )
  {
    ctxt_.ibcast_send(1,1,&sp.atomic_number_,1);
    ctxt_.dbcast_send(1,1,&sp.mass_,1);
    ctxt_.ibcast_send(1,1,&sp.zval_,1);
    ctxt_.ibcast_send(1,1,&sp.lmax_,1);
    ctxt_.ibcast_send(1,1,&sp.llocal_,1);
    ctxt_.ibcast_send(1,1,&sp.nquad_,1);
    ctxt_.dbcast_send(1,1,&sp.rquad_,1);
    ctxt_.dbcast_send(1,1,&sp.deltar_,1);
  }
  else
  {
    ctxt_.ibcast_recv(1,1,&sp.atomic_number_,1,0,0);
    ctxt_.dbcast_recv(1,1,&sp.mass_,1,0,0);
    ctxt_.ibcast_recv(1,1,&sp.zval_,1,0,0);
    ctxt_.ibcast_recv(1,1,&sp.lmax_,1,0,0);
    ctxt_.ibcast_recv(1,1,&sp.llocal_,1,0,0);
    ctxt_.ibcast_recv(1,1,&sp.nquad_,1,0,0);
    ctxt_.dbcast_recv(1,1,&sp.rquad_,1,0,0);
    ctxt_.dbcast_recv(1,1,&sp.deltar_,1,0,0);
    
    sp.vps_.resize(sp.lmax_+1);
    sp.phi_.resize(sp.lmax_+1);
  }
  
  ctxt_.string_bcast(sp.symbol_,0);
  ctxt_.string_bcast(sp.description_,0);
  ctxt_.string_bcast(sp.uri_,0);
  
  for ( int l = 0; l <= sp.lmax_; l++ )
  {
    int np_vps;
    if ( ctxt_.onpe0() )
    {
      np_vps = sp.vps_[l].size();
      ctxt_.ibcast_send(1,1,&np_vps,1);
      ctxt_.dbcast_send(np_vps,1,&sp.vps_[l][0],np_vps);
    }
    else
    {
      ctxt_.ibcast_recv(1,1,&np_vps,1,0,0);
      sp.vps_[l].resize(np_vps);
      ctxt_.dbcast_recv(np_vps,1,&sp.vps_[l][0],np_vps,0,0);
    }
  }

  // broadcast atomic orbitals
  for ( int l = 0; l <= sp.lmax_; l++ )
  {
    int np_phi;
    if ( ctxt_.onpe0() )
    {
      np_phi = sp.phi_[l].size();
      ctxt_.ibcast_send(1,1,&np_phi,1);
      ctxt_.dbcast_send(np_phi,1,&sp.phi_[l][0],np_phi);
    }
    else
    {
      ctxt_.ibcast_recv(1,1,&np_phi,1,0,0);
      sp.phi_[l].resize(np_phi);
      ctxt_.dbcast_recv(np_phi,1,&sp.phi_[l][0],np_phi,0,0);
    }
  }
}
