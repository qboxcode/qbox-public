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
// B3LYPFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef B3LYPFUNCTIONAL_H
#define B3LYPFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>

class B3LYPFunctional : public XCFunctional
{
  B3LYPFunctional();

  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn,
                 _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];

  void excb3lyp(double rho, double grad,
    double *exc, double *vxc1, double *vxc2);

  void excb3lyp_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *exc_up, double *exc_dn, double *vxc1_up, double *vxc1_dn,
    double *vxc2_upup, double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup);

  public:

  B3LYPFunctional(const std::vector<std::vector<double> > &rhoe);

  bool isGGA() const { return true; };
  std::string name() const { return "B3LYP"; };
  void setxc(void);
};
#endif
