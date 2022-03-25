//////////////////////////////////////////////////////////////////////////////// //
// Copyright (c) 2018 The Regents of the University of California
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
// Function3dHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Function3d.h"
#include "Function3dHandler.h"
#include "Base64Transcoder.h"
#include "StrX.h"
#include "Timer.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
using namespace xercesc;
#include <iostream>
#include <sstream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Function3dHandler::Function3dHandler(Function3d& f) : f_(f) {}

////////////////////////////////////////////////////////////////////////////////
Function3dHandler::~Function3dHandler(void) {}

////////////////////////////////////////////////////////////////////////////////
void Function3dHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  string locname = StrX(localname).localForm();
  if ( locname == "function3d" )
  {
    unsigned int len = attributes.getLength();
    for ( unsigned int index = 0; index < len; index++ )
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "name" )
      {
        f_.name = StrX(attributes.getValue(index)).localForm();
      }
    }
  }
  else if ( locname == "domain" )
  {
    unsigned int len = attributes.getLength();
    for ( unsigned int index = 0; index < len; index++ )
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "a" )
      {
        istringstream stst(StrX(attributes.getValue(index)).localForm());
        stst >> f_.a;
      }
      else if ( attrname == "b" )
      {
        istringstream stst(StrX(attributes.getValue(index)).localForm());
        stst >> f_.b;
      }
      else if ( attrname == "c" )
      {
        istringstream stst(StrX(attributes.getValue(index)).localForm());
        stst >> f_.c;
      }
    }
  }
  else if ( locname == "grid" )
  {
    unsigned int len = attributes.getLength();
    for ( unsigned int index = 0; index < len; index++ )
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "nx" )
      {
        f_.nx = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "ny" )
      {
        f_.ny = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nz" )
      {
        f_.nz = atoi(StrX(attributes.getValue(index)).localForm());
      }
    }
  }
  else if ( locname == "grid_function" )
  {
    x0_ = 0;
    y0_ = 0;
    z0_ = 0;
    unsigned int len = attributes.getLength();
    for ( unsigned int index = 0; index < len; index++ )
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "type" )
      {
        string type = StrX(attributes.getValue(index)).localForm();
        if ( type != "double" )
        {
          cerr << "Function3dHandler: grid_function type must be double"
               << endl;
          assert(!"Function3dHandler: incorrect grid_function type");
        }
      }
      else if ( attrname == "nx" )
      {
        fnx_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "ny" )
      {
        fny_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nz" )
      {
        fnz_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "x0" )
      {
        x0_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "y0" )
      {
        y0_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "z0" )
      {
        z0_ = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "encoding" )
      {
        string type = StrX(attributes.getValue(index)).localForm();
        if ( type != "base64" )
        {
          cerr << "Function3dHandler: encoding must be base64"
               << endl;
          assert(!"Function3dHandler: incorrect encoding");
        }
      }
    }
    buf_ = "";
  }
}

////////////////////////////////////////////////////////////////////////////////
#ifndef XERCES_VERSION_MAJOR
#error "XERCES_VERSION_MAJOR not defined"
#endif
#if XERCES_VERSION_MAJOR < 3
#error Xerces-C version should be at least 3
#endif
void Function3dHandler::characters(const XMLCh* const chars,
  const XMLSize_t length)
{
#if TIMING
  Timer tm;
  tm.start();
#endif

  char *str = XMLString::transcode(chars);
  buf_ += str;
  XMLString::release(&str);

#if TIMING
  tm.stop();
  cout << " Function3dHandler::characters: time: " << tm.real() << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Function3dHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname)
{
  string locname = StrX(localname).localForm();
  if ( locname == "grid_function" )
  {
    if ( (fnx_ > f_.nx) || (fny_ > f_.ny) || (fnz_ > f_.nz) )
    {
      assert(!"Function3dHandler:: fragment size > grid size");
    }
    if ( (fnx_ != f_.nx) || (fny_ != f_.ny) || (fnz_ != f_.nz) )
    {
      assert(!"Function3dHandler:: fragment processing not implemented");
    }
    if ( x0_ != 0 || y0_ != 0 || z0_ != 0 )
    {
      assert(!"Function3dHandler:: fragment offset not implemented");
    }
    // process base64 data in buf
    f_.val.resize(f_.nx*f_.ny*f_.nz);
    Timer tm;
    tm.start();
    Base64Transcoder xcdr;
    size_t nbytes = xcdr.decode(buf_.size(),buf_.data(),(byte_t*)&f_.val[0]);
    assert(nbytes==f_.val.size()*sizeof(double));
    buf_.clear();
    #if PLT_BIG_ENDIAN
    xcdr.byteswap_double(f_.val.size(),&f_.val[0]);
    #endif
  }
}
