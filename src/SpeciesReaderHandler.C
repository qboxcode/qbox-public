////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReaderHandler.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReaderHandler.C,v 1.1 2003-03-27 22:05:59 fgygi Exp $

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include "SpeciesReaderHandler.h"
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SpeciesReaderHandler::SpeciesReaderHandler(Species& sp) : 
  sp_(sp) {}

////////////////////////////////////////////////////////////////////////////////
SpeciesReaderHandler::~SpeciesReaderHandler() {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::startDocument() {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{  
  st = "";

  unsigned int len = attributes.getLength();
  for (unsigned int index = 0; index < len; index++)
  {
    string attrname(XMLString::transcode(attributes.getLocalName(index)));
    if ( attrname == "l")
    {
      current_l = atoi(StrX(attributes.getValue(index)).localForm());
    }
    else if ( attrname == "size" )
    {
      current_size = atoi(StrX(attributes.getValue(index)).localForm());
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::characters(const XMLCh* const chars,
  const unsigned int length)
{
    st += string(XMLString::transcode(chars),length);
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname)
{
  istringstream stst(st);
  string locname(XMLString::transcode(localname));
  if ( locname == "description")
  {
    sp_.description_ = st;
  }
  else if ( locname == "atomic_number" )
  {
    stst >> sp_.atomic_number_;
  }
  else if ( locname == "mass" )
  {
    stst >> sp_.mass_;
  }
  else if ( locname == "symbol" )
  {
    stst >> skipws >> sp_.symbol_;
  }
  else if ( locname == "valence_charge" )
  {
    stst >> sp_.zval_;
  }
  else if ( locname == "lmax" )
  {
    stst >> sp_.lmax_;
  }
  else if ( locname == "llocal" )
  {
    stst >> sp_.llocal_;
  }
  else if ( locname == "nquad" )
  {
    stst >> sp_.nquad_;
  }
  else if ( locname == "rquad" )
  {
    stst >> sp_.rquad_;
  }
  else if ( locname == "mesh_spacing" )
  {
    stst >> sp_.deltar_;
  }
  else if ( locname == "radial_potential" )
  {
    if ( current_l+1  > sp_.vps_.size() )
    {
      sp_.vps_.resize(current_l+1);
      sp_.phi_.resize(current_l+1);
    }
    sp_.vps_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.vps_[current_l][i];
  }
  else if ( locname == "radial_function" )
  {
    sp_.phi_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.phi_[current_l][i];
  }
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::endDocument() {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::ignorableWhitespace( const   XMLCh* const chars,
  const  unsigned int length) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::processingInstruction(const  XMLCh* const target,
  const XMLCh* const data) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::error(const SAXParseException& e)
{
    cout << "\nError at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
         
    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::fatalError(const SAXParseException& e)
{
    cout << "\nFatal Error at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::warning(const SAXParseException& e)
{
    cout << "\nWarning at file " << StrX(e.getSystemId())
		 << ", line " << e.getLineNumber()
		 << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::unparsedEntityDecl(const     XMLCh* const name,
  const   XMLCh* const publicId, const   XMLCh* const systemId,
  const   XMLCh* const notationName)
{
    // Not used at this time
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReaderHandler::notationDecl(const   XMLCh* const name,
  const XMLCh* const publicId, const XMLCh* const systemId)
{
    // Not used at this time
}
