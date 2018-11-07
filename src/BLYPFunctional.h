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
// BLYPFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
#ifndef BLYPFUNCTIONAL_H
#define BLYPFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>

class BLYPFunctional : public XCFunctional
{
  BLYPFunctional();

  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn,
                 _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];

  public:

  BLYPFunctional(const std::vector<std::vector<double> > &rhoe);

  static void exb88(double rho, double grad,
    double *ex, double *vx1, double *vx2);
  static void eclyp(double rho, double grad,
    double *ec, double *vc1, double *vc2);

  static void exb88_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *ex_up, double *ex_dn, double *vx1_up, double *vx1_dn,
    double *vx2_upup, double *vx2_dndn, double *vx2_updn, double *vx2_dnup);
  static void eclyp_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *ec_up, double *ec_dn, double *vc1_up, double *vc1_dn,
    double *vc2_upup, double *vc2_dndn, double *vc2_updn, double *vc2_dnup);

  bool isGGA() const { return true; };
  std::string name() const { return "BLYP"; };
  void setxc(void);
};
#endif
