////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionHandler.h,v 1.9 2007-10-19 16:24:05 fgygi Exp $

#ifndef WavefunctionHANDLER_H
#define WavefunctionHANDLER_H

#include "StructureHandler.h"
#include "UnitCell.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class FourierTransform;

class WavefunctionHandler : public StructureHandler
{
  private:

  Wavefunction& wf_;
  DoubleMatrix& gfdata_;
  UnitCell uc;
  UnitCell ruc;
  double ecut;
  std::vector<double> dmat;
  int nx,ny,nz;
  int current_gf_nx,current_gf_ny,current_gf_nz;
  std::string current_gf_encoding;
  int current_ispin,current_ikp,current_n,current_igf;
  int read_from_gfdata;
  FourierTransform* ft;
  std::vector<std::complex<double> > wftmp;

  void byteswap_double(size_t n, double* x);

  public:

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

  WavefunctionHandler(Wavefunction& wf, DoubleMatrix& gfdata);
  ~WavefunctionHandler();
};
#endif
