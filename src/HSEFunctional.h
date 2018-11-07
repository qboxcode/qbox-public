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
// HSEFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Implements the HSE hybrid functional
// Heyd et al.,   J. Chem. Phys. 118, 8207 (2003)
// Heyd et al.,   J. Chem. Phys. 120, 7274 (2004)
// Krukau et al., J. Chem. Phys. 125, 224106 (2006)
//
////////////////////////////////////////////////////////////////////////////////

#ifndef HSEFUNCTIONAL_H
#define HSEFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>

class HSEFunctional : public XCFunctional
{
  const double x_coeff_, c_coeff_;

  // screening parameter of the HSE functional
  const double omega; // == 0.11

  // vectors common to all GGA exchange functionals
  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn, _vxc2, _vxc2_upup, _vxc2_updn,
    _vxc2_dnup, _vxc2_dndn;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];

  public:

  // constructor
  HSEFunctional(const std::vector<std::vector<double> > &rhoe);

  // HSE's local part is derived from PBE
  bool isGGA() const { return true; }

  // return the name of the functional
  std::string name() const { return "HSE"; }

  void setxc(void);
};
#endif
