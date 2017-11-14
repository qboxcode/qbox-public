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
// SCANFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SCANFUNCTIONAL_H
#define SCANFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>

class SCANFunctional : public XCFunctional
{
  double x_coeff_, c_coeff_;
  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn,
                 _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn,
                 _vxc3;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];
  std::vector<double> _tau;

  void gPW92(double alpha, double beta0, double beta1, double beta2,
    double beta3, double beta4, double rtrs, double *gg, double *dgdrs);

  void excSCAN(double rho, double grad, double tau, double *exc,
    double *vxc1, double *vxc2, double *vxc3);

  void excSCAN_sp(double rho_up, double rho_dn,
    double grad_up, double grad_dn, double grad,
    double *exc_up, double *exc_dn,
    double *vxc1_up, double *vxc1_dn, double *vxc2_upup, double *vxc2_dndn,
    double *vxc2_updn, double *vxc2_dnup);

  public:

  // constructor with variable coefficients for exchange and correlation
  // with default values 1.0
  SCANFunctional(const std::vector<std::vector<double> > &rhoe,
                double x_coeff=1.0, double c_coeff=1.0);

  bool isGGA() const { return true; };
  bool isMeta() const { return true; };
  std::string name() const { return "SCAN"; };
  void setxc(void);
};
#endif
