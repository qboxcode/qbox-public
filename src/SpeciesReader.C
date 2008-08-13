////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 2 of 
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.C,v 1.10 2008-08-13 06:39:43 fgygi Exp $

#include "Species.h"
#include "SpeciesReader.h"
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#if USE_XERCES
#include "StructuredDocumentHandler.h"
#include "SpeciesHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
using namespace xercesc;
#else
#include <sstream>
#include <cstdio>
#include <sys/stat.h>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SpeciesReader::SpeciesReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::readSpecies (Species& sp, const string uri)
{
#if USE_XERCES
  if ( ctxt_.onpe0() )
  {
    const char* encodingName = "UTF-8";
    SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
    //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
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
      cout << "  Species::readSpecies: Error during XML initialization :\n"
           << StrX(toCatch.getMessage()) << endl;
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

#else

  // USE_XERCES is not defined

  // Simplified parsing not using Xerces
  // parse file uri

  if ( ctxt_.onpe0() )
  {
    struct stat statbuf;
    bool found_file = !stat(uri.c_str(),&statbuf);
    assert(found_file);
      cout << "  SpeciesReader opening file "
           << uri << " size: "
           << statbuf.st_size << endl;

    FILE* infile;
    infile = fopen(uri.c_str(),"r");
    if ( !infile )
    {
      cout << "  SpeciesReader::readSpecies could not open file "
           << uri << " for reading" << endl;
      return;
    }
    off_t sz = statbuf.st_size;
    string buf;
    buf.resize(sz);
    fread(&buf[0],sizeof(char),sz,infile);

    string::size_type pos = 0;

    string tag, start_tag, end_tag;
    string::size_type start, end, len;

    tag = "description";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;

    sp.description_ = buf.substr(start,len);
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.description_ << endl;

    tag = "symbol";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.symbol_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.symbol_ << endl;

    tag = "atomic_number";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.atomic_number_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.atomic_number_ << endl;

    tag = "mass";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.mass_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.mass_ << endl;

    tag = "valence_charge";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.zval_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.zval_ << endl;

    tag = "lmax";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.lmax_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.lmax_ << endl;

    tag = "llocal";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.llocal_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.llocal_ << endl;

    tag = "nquad";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.nquad_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.nquad_ << endl;

    tag = "rquad";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.rquad_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.rquad_ << endl;

    tag = "mesh_spacing";
    start_tag = string("<") + tag + string(">");
    end_tag = string("</") + tag + string(">");
    start = buf.find(start_tag,pos);
    assert(start != string::npos );
    start = buf.find(">",start)+1;
    end = buf.find(end_tag);
    pos = buf.find(">",end)+1;
    len = end - start;
    {
      istringstream stst(buf.substr(start,len));
      stst >> sp.deltar_;
    }
    cout << "  SpeciesReader::readSpecies: read " << tag << " "
         << sp.deltar_ << endl;

    for ( int l = 0; l < sp.lmax_ + 1; l++ )
    {
      // read projector
      int size;
      tag = "projector";
      start_tag = string("<") + tag;
      start = buf.find(start_tag,pos);
      assert(start != string::npos );

      pos = buf.find("l=",pos)+3;
      start = pos;
      end = buf.find("\"",start);
      len = end - start;

      {
        istringstream stst(buf.substr(start,len));
        int lread;
        stst >> lread;
        //cout << " lread=" << lread << endl;
        //cout << " l=" << l << endl;
        assert(l==lread);
      }

      pos = buf.find("size=",pos)+6;
      start = pos;
      end = buf.find("\"",start);
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        stst >> size;
        // cout << " size=" << size << endl;
      }
      // read radial potential

      sp.vps_.resize(sp.vps_.size()+1);
      sp.vps_[l].resize(size);

      tag = "radial_potential";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      assert(start != string::npos );
      start = buf.find(">",start)+1;
      end = buf.find(end_tag,start);
      pos = buf.find(">",end)+1;
      len = end - start;
      {
        istringstream stst(buf.substr(start,len));
        for ( int i = 0; i < size; i++ )
        {
          stst >> sp.vps_[l][i];
          //cout << sp.vps_[l][i] << endl;
        }
      }
      cout << "  SpeciesReader::readSpecies: read " << tag << " l="
           << l << " size=" << size << endl;

      sp.phi_.resize(sp.phi_.size()+1);
      sp.phi_[l].resize(size);

      // read radial function only if the radial_function tag was found

      tag = "radial_function";
      start_tag = string("<") + tag + string(">");
      end_tag = string("</") + tag + string(">");
      start = buf.find(start_tag,pos);
      if ( l != sp.llocal_ )
      {
        // if l is not the local potential, there must be a radial function
        assert(start != string::npos );
      }

      if ( start != string::npos )
      {
        start = buf.find(">",start)+1;
        end = buf.find(end_tag,start);
        pos = buf.find(">",end)+1;
        len = end - start;
        {
          istringstream stst(buf.substr(start,len));
          for ( int i = 0; i < size; i++ )
          {
            stst >> sp.phi_[l][i];
            //cout << sp.phi_[l][i] << endl;
          }
        }
        cout << "  SpeciesReader::readSpecies: read " << tag << " l="
             << l << " size=" << size << endl;
      }
    }

    sp.uri_ = uri;
  }
#endif
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
