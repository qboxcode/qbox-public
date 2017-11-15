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
// SampleHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLEHANDLER_H
#define SAMPLEHANDLER_H

#include "StructureHandler.h"
#include <string>
#include <vector>
class DoubleMatrix;
class Sample;
class Wavefunction;

class SampleHandler : public StructureHandler
{
  private:

  Sample& s_;
  DoubleMatrix& gfdata_;
  int current_gfdata_pos;

  public:

  bool read_wf,read_wfv;

  // Start of the root element in the structure being handled
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of the root element in the structure being handled
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  SampleHandler(Sample& s, DoubleMatrix& gfdata);
  ~SampleHandler();
};
#endif
