////////////////////////////////////////////////////////////////////////////////
//
// XCPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XCPotential.h,v 1.5 2007-10-19 16:24:05 fgygi Exp $

#ifndef XCPOTENTIAL_H
#define XCPOTENTIAL_H

#include "ChargeDensity.h"
#include "LDAFunctional.h"
#include "PBEFunctional.h"
#include "BLYPFunctional.h"
#include <string>
#include <vector>
#include <valarray>
#include <complex>

class Basis;
class FourierTransform;

class XCPotential
{
  private:

  const Context& ctxt_;
  const ChargeDensity& cd_;
  XCFunctional* xcf_;

  std::vector<std::vector<double> > vxctmp;          // vxctmp[ispin][ir]
  std::vector<std::complex<double> > tmpr;           // tmpr[ir]
  std::vector<std::complex<double> > tmp1, tmp2;     // tmp1[ig], tmp2[ig]

  double exc_, dxc_, dxc0_, dxc1_, dxc2_;
  int nspin_;
  int ngloc_;
  int np012loc_;

  FourierTransform& vft_;
  Basis& vbasis_;

  public:

  const XCFunctional* xcf() { return xcf_; }
  XCPotential(const ChargeDensity& cd, const std::string functional_name);
  ~XCPotential();
  void update(std::vector<std::vector<double> >& vr);
  void compute_stress(std::valarray<double>& sigma_exc);
  double exc(void) { return exc_; }
};

class XCPotentialException
{
  public:
  std::string msg;
  XCPotentialException(std::string s) : msg(s) {}
};
#endif
