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
// RSHFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Range-separated hybrid functional (RSH)
// J.H.Skone et al. Phys. Rev. B93, 235106 (2016)
// RSH is defined by alpha_RSH, beta_RSH, mu_RSH
// sigma = alpha_RSH * rho(r,r') * erf(r-r')/(r-r') +
//         beta_RSH * rho(r,r') * erfc(r-r')/(r-r') +
//         (1 - alpha_RSH) * Vx_LR(r,mu_RSH) +
//         (1 - beta_RSH) * Vx_SR(r,mu_RSH)
// The HSE functional is obtained using alpha_RSH=0, beta_RSH=0.25, mu_RSH=0.11
// Heyd et al.,   J. Chem. Phys. 118, 8207 (2003)
// Heyd et al.,   J. Chem. Phys. 120, 7274 (2004)
// Krukau et al., J. Chem. Phys. 125, 224106 (2006)
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RSHFUNCTIONAL_H
#define RSHFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>

class RSHFunctional : public XCFunctional
{
  const double x_coeff_, c_coeff_;
  const double alpha_RSH_, beta_RSH_, mu_RSH_;
  const double omega; // == mu_RSH

  // vectors common to all GGA exchange functionals
  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn, _vxc2, _vxc2_upup, _vxc2_updn,
    _vxc2_dnup, _vxc2_dndn;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];

  void RSH_exchange(const double rho, const double grad,
  const double a_ex, const double w, double *ex, double *vx1, double *vx2);

  void gcor2(double a, double a1, double b1, double b2,
  double b3, double b4, double rtrs, double *gg, double *ggrs);

  void PBE_correlation(const double rho, const double grad,
  double *ec, double *vc1, double *vc2);

  void PBE_correlation_sp(const double rho_up, const double rho_dn,
  const double grad_up, const double grad_dn, const double grad, double *ec,
  double *vc1_up, double *vc1_dn, double *vc2);

  void approximateIntegral(const double omega_kF, const double Hs2,
  const double D_term, const double dHs2_ds, double *appInt,
  double *dAppInt_ds, double *dAppInt_dkF);

  void RSH_enhance(const double s_inp, const double kF,
  const double w, double *fx, double *dfx_ds, double* dfx_dkf);

  public:

  // constructor
  RSHFunctional(const std::vector<std::vector<double> > &rhoe,
                double alpha_RSH, double beta_RSH, double mu_RSH);

  bool isGGA() const { return true; }

  // return the name of the functional
  std::string name() const { return "RSH"; }

  void setxc(void);
};

#endif
