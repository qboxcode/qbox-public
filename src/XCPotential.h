////////////////////////////////////////////////////////////////////////////////
//
// XCPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XCPotential.h,v 1.1 2003-04-10 19:16:52 fgygi Exp $

#ifndef XCPOTENTIAL_H
#define XCPOTENTIAL_H

#include "ChargeDensity.h"
#include "LDAFunctional.h"
#include "PBEFunctional.h"
#include "BLYPFunctional.h"
#include <vector>
#include <complex>
using namespace std;

class Basis;
class FourierTransform;

class XCPotential
{
  private:
  
  const Context& ctxt_;  
  ChargeDensity& cd_;
  XCFunctional* xcf_;
  
  vector<vector<double> > vxctmp;          // vxctmp[ispin][ir]
  vector<complex<double> > tmpr;           // tmpr[ir]
  vector<complex<double> > tmp1, tmp2;     // tmp1[ig], tmp2[ig]
  
  double exc_, dxc_, dxc0_, dxc1_, dxc2_;
  int nspin_;
  int ngloc_;
  int np012loc_;
  
  FourierTransform& vft_;
  Basis& vbasis_;

  public:

  const XCFunctional* xcf() { return xcf_; }
  XCPotential(ChargeDensity& cd, const string functional_name);
  ~XCPotential();
  void update(vector<vector<double> >& vr, bool compute_stress);
  double exc(void) { return exc_; }
  double dxc(void) { return dxc_; }
  double dxc0(void) { return dxc0_; }
  double dxc1(void) { return dxc1_; }
  double dxc2(void) { return dxc2_; }
};

class XCPotentialException
{
  public:
  string msg;
  XCPotentialException(string s) : msg(s) {}
};
#endif
