////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2015 The Regents of the University of California
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
// VWNFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VWNFUNCTIONAL_H
#define VWNFUNCTIONAL_H

#include <vector>
#include <cassert>
#include "XCFunctional.h"

class VWNFunctional : public XCFunctional
{
  static void alpha_c(const double rh, double &a, double &da);
  static void x_unpolarized(const double rh, double &ex, double &vx);
  static void c_unpolarized(const double rh, double &ec, double &vc);
  static void x_polarized(const double rh, double &ex, double &vx);
  static void c_polarized(const double rh, double &ec, double &vc);

  std::vector<double> _exc;
  std::vector<std::vector<double> > _vxc;

  VWNFunctional();

  public:

  VWNFunctional(const std::vector<std::vector<double> > &rhoe)
  {
    _nspin = rhoe.size();
    if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
    _np = rhoe[0].size();
    _exc.resize(_np);
    _vxc.resize(_nspin);
    for ( int i = 0; i < _nspin; i++ )
    {
      _vxc[i].resize(_np);
    }

    if ( _nspin == 1 )
    {
      rho = &rhoe[0][0];
      exc = &_exc[0];
      vxc1 = &_vxc[0][0];
    }
    else
    {
      rho_up = &rhoe[0][0];
      rho_dn = &rhoe[1][0];
      exc = &_exc[0];
      vxc1_up = &_vxc[0][0];
      vxc1_dn = &_vxc[1][0];
    }
  };

  static void exvwn(const double rh, double &ex, double &vx);
  static void ecvwn(const double rh, double &ec, double &vc);

  static void exvwn_sp(double roe_up, double roe_dn,
    double &ex, double &vx_up, double &vx_dn);
  static void ecvwn_sp(double roe_up, double roe_dn,
    double &ec, double &vc_up, double &vc_dn);

  std::string name() const { return "VWN"; };
  void setxc(void);
};
#endif
