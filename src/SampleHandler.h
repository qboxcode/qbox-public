////////////////////////////////////////////////////////////////////////////////
//
// SampleHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleHandler.h,v 1.6 2007-10-19 16:24:04 fgygi Exp $

#ifndef SAMPLEHANDLER_H
#define SAMPLEHANDLER_H

#include "StructureHandler.h"
#include <string>
class DoubleMatrix;
class Sample;
class Wavefunction;

class SampleHandler : public StructureHandler
{
  private:

  Sample& s_;
  DoubleMatrix& gfdata_;
  Wavefunction& wfvtmp_;

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

  SampleHandler(Sample& s, DoubleMatrix& gfdata, Wavefunction& wfvtmp);
  ~SampleHandler();
};
#endif
