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
// SampleHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SampleHandler.h"
#include "Sample.h"
#include "AtomSetHandler.h"
#include "WavefunctionHandler.h"
#include "StrX.h"
using namespace xercesc;
#include <iostream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleHandler::SampleHandler(Sample& s, DoubleMatrix& gfdata) :
  s_(s), gfdata_(gfdata), current_gfdata_pos(0) {}

////////////////////////////////////////////////////////////////////////////////
SampleHandler::~SampleHandler() {}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " SampleHandler::startElement " << StrX(qname) << endl;
}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  // istringstream stst(st);
  // string locname = StrX(localname).localForm();
  // cout << " SampleHandler::endElement " << locname << endl;
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* SampleHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " SampleHandler::startSubHandler " << StrX(qname) << endl;

  string qnm = StrX(qname).localForm();
  if ( qnm == "atomset" )
  {
    return new AtomSetHandler(s_.atoms);
  }
  else if ( qnm == "wavefunction" )
  {
    read_wf = true;
    return new WavefunctionHandler(s_.wf,gfdata_,current_gfdata_pos);
  }
  else if ( qnm == "wavefunction_velocity" )
  {
    read_wfv = true;
    assert(read_wf);
    s_.wfv = new Wavefunction(s_.wf);
    // reset wfv to have only k=0
    s_.wfv->reset();
    return new WavefunctionHandler(*s_.wfv,gfdata_,current_gfdata_pos);
  }
  else
  {
    return 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler)
{
  // cout << " SampleHandler::endSubHandler " << StrX(qname) << endl;
  delete subHandler;
}
