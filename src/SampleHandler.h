////////////////////////////////////////////////////////////////////////////////
//
// SampleHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleHandler.h,v 1.4 2003-09-23 19:03:21 fgygi Exp $

#ifndef SAMPLEHANDLER_H
#define SAMPLEHANDLER_H

#include "StructureHandler.h"
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
      const XMLCh* const qname, string& content);
  
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
