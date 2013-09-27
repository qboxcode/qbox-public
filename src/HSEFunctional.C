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
// HSEFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
//
// Implements the HSE hybrid functional
// Heyd et al.,   J. Chem. Phys. 118, 8207 (2003)
// Heyd et al.,   J. Chem. Phys. 120, 7274 (2004)
// Krukau et al., J. Chem. Phys. 125, 224106 (2006)
//
////////////////////////////////////////////////////////////////////////////////
//
// Author: Martin Schlipf (2013)
// Contact: martin.schlipf@gmail.com
//

#include "HSEFunctional.h"
#include <cassert>

using namespace std;

// constructor
HSEFunctional::HSEFunctional(const vector<vector<double> > &rhoe) :
  x_coeff_(0.75), c_coeff_(1.0)
{
  // nonmagnetic or magnetic
  _nspin = rhoe.size();
  // in magnetic calculations both density arrays must have the same size
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());

  // allocate arrays
  _np = rhoe[0].size();
  if ( _nspin == 1 )
  {
    // nonmagnetic arrays used
    _exc.resize(_np);
    _vxc1.resize(_np);
    _vxc2.resize(_np);
    _grad_rho[0].resize(_np);
    _grad_rho[1].resize(_np);
    _grad_rho[2].resize(_np);
    rho = &rhoe[0][0];
    grad_rho[0] = &_grad_rho[0][0];
    grad_rho[1] = &_grad_rho[1][0];
    grad_rho[2] = &_grad_rho[2][0];
    exc = &_exc[0];
    vxc1 = &_vxc1[0];
    vxc2 = &_vxc2[0];
  }
  else
  {
    // magnetic arrays used
    _exc_up.resize(_np);
    _exc_dn.resize(_np);
    _vxc1_up.resize(_np);
    _vxc1_dn.resize(_np);
    _vxc2_upup.resize(_np);
    _vxc2_updn.resize(_np);
    _vxc2_dnup.resize(_np);
    _vxc2_dndn.resize(_np);
    _grad_rho_up[0].resize(_np);
    _grad_rho_up[1].resize(_np);
    _grad_rho_up[2].resize(_np);
    _grad_rho_dn[0].resize(_np);
    _grad_rho_dn[1].resize(_np);
    _grad_rho_dn[2].resize(_np);

    rho_up = &rhoe[0][0];
    rho_dn = &rhoe[1][0];
    grad_rho_up[0] = &_grad_rho_up[0][0];
    grad_rho_up[1] = &_grad_rho_up[1][0];
    grad_rho_up[2] = &_grad_rho_up[2][0];
    grad_rho_dn[0] = &_grad_rho_dn[0][0];
    grad_rho_dn[1] = &_grad_rho_dn[1][0];
    grad_rho_dn[2] = &_grad_rho_dn[2][0];
    exc_up = &_exc_up[0];
    exc_dn = &_exc_dn[0];
    vxc1_up = &_vxc1_up[0];
    vxc1_dn = &_vxc1_dn[0];
    vxc2_upup = &_vxc2_upup[0];
    vxc2_updn = &_vxc2_updn[0];
    vxc2_dnup = &_vxc2_dnup[0];
    vxc2_dndn = &_vxc2_dndn[0];
  }

}

void HSEFunctional::setxc(void)
{
  // dummy
}
