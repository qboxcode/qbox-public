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
// BHandHLYPFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include "BHandHLYPFunctional.h"
#include "BLYPFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BHandHLYPFunctional::BHandHLYPFunctional(const vector<vector<double> > &rhoe)
{
  _nspin = rhoe.size();
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
  _np = rhoe[0].size();

  if ( _nspin == 1 )
  {
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

////////////////////////////////////////////////////////////////////////////////
void BHandHLYPFunctional::setxc(void)
{
  if ( _np == 0 ) return;

  if ( _nspin == 1 )
  {
    assert( rho != 0 );
    assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
    assert( exc != 0 );
    assert( vxc1 != 0 );
    assert( vxc2 != 0 );
    for ( int i = 0; i < _np; i++ )
    {
      double grad = sqrt(grad_rho[0][i]*grad_rho[0][i] +
                         grad_rho[1][i]*grad_rho[1][i] +
                         grad_rho[2][i]*grad_rho[2][i] );
      excbhandhlyp(rho[i],grad,&exc[i],&vxc1[i],&vxc2[i]);
    }
  }
  else
  {
    assert( rho_up != 0 );
    assert( rho_dn != 0 );
    assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
    assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
    assert( exc_up != 0 );
    assert( exc_dn != 0 );
    assert( vxc1_up != 0 );
    assert( vxc1_dn != 0 );
    assert( vxc2_upup != 0 );
    assert( vxc2_updn != 0 );
    assert( vxc2_dnup != 0 );
    assert( vxc2_dndn != 0 );

    for ( int i = 0; i < _np; i++ )
    {
      double grx_up = grad_rho_up[0][i];
      double gry_up = grad_rho_up[1][i];
      double grz_up = grad_rho_up[2][i];
      double grx_dn = grad_rho_dn[0][i];
      double gry_dn = grad_rho_dn[1][i];
      double grz_dn = grad_rho_dn[2][i];
      double grad_up2 = grx_up*grx_up + gry_up*gry_up + grz_up*grz_up;
      double grad_dn2 = grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn;
      double grad_up_grad_dn = grx_up*grx_dn + gry_up*gry_dn + grz_up*grz_dn;

      excbhandhlyp(rho_up[i],rho_dn[i],grad_up2,grad_dn2,grad_up_grad_dn,
                  &exc_up[i],&exc_dn[i],&vxc1_up[i],&vxc1_dn[i],
                  &vxc2_upup[i],&vxc2_dndn[i],
                  &vxc2_updn[i],&vxc2_dnup[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void BHandHLYPFunctional::excbhandhlyp(double rho, double grad,
  double *exc, double *vxc1, double *vxc2)
{
  // BHandHLYP unpolarized
  // Coefficients of the BHandHLYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // B3LYP: 0.2 ExHF + 0.8 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  // BHandHLYP: 0.5 ExHF + 0.5 ExB88 + 1.0 EcLYP
  const double xb88_coeff = 0.5; // Becke88 exchange
  const double clyp_coeff = 1.0 ;

  double ex_b88,vx1_b88,vx2_b88;
  double ec_lyp,vc1_lyp,vc2_lyp;

  BLYPFunctional::exb88(rho,grad,&ex_b88,&vx1_b88,&vx2_b88);
  BLYPFunctional::eclyp(rho,grad,&ec_lyp,&vc1_lyp,&vc2_lyp);

  *exc = xb88_coeff  * ex_b88 +
         clyp_coeff  * ec_lyp;

  *vxc1 = xb88_coeff * vx1_b88 +
          clyp_coeff * vc1_lyp;

  *vxc2 = xb88_coeff * vx2_b88 +
          clyp_coeff * vc2_lyp;
}

////////////////////////////////////////////////////////////////////////////////
void BHandHLYPFunctional::excbhandhlyp(double rho_up, double rho_dn,
  double grad_up2, double grad_dn2, double grad_up_grad_dn,
  double *exc_up, double *exc_dn, double *vxc1_up, double *vxc1_dn,
  double *vxc2_upup, double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup)
{
  // BHandHLYP spin-polarized
  // Coefficients of the BHandHLYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // B3LYP: 0.2 ExHF + 0.8 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  // BHandHLYP: 0.5 ExHF + 0.5 ExB88 + 1.0 EcLYP
  const double xb88_coeff = 0.5; // Becke88 exchange
  const double clyp_coeff = 1.0;

  double ex_b88_up,ex_b88_dn,vx1_b88_up,vx1_b88_dn,
         vx2_b88_upup,vx2_b88_dndn,vx2_b88_updn,vx2_b88_dnup;
  double ec_lyp_up,ec_lyp_dn,vc1_lyp_up,vc1_lyp_dn,
         vc2_lyp_upup,vc2_lyp_dndn,vc2_lyp_updn,vc2_lyp_dnup;

  BLYPFunctional::exb88_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ex_b88_up,&ex_b88_dn,&vx1_b88_up,&vx1_b88_dn,
    &vx2_b88_upup,&vx2_b88_dndn,&vx2_b88_updn,&vx2_b88_dnup);
  BLYPFunctional::eclyp_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ec_lyp_up,&ec_lyp_dn,&vc1_lyp_up,&vc1_lyp_dn,
    &vc2_lyp_upup,&vc2_lyp_dndn,&vc2_lyp_updn,&vc2_lyp_dnup);

  *exc_up = xb88_coeff  * ex_b88_up +
            clyp_coeff  * ec_lyp_up;

  *exc_dn = xb88_coeff  * ex_b88_dn +
            clyp_coeff  * ec_lyp_dn;

  *vxc1_up = xb88_coeff * vx1_b88_up +
             clyp_coeff * vc1_lyp_up;

  *vxc1_dn = xb88_coeff * vx1_b88_dn +
             clyp_coeff * vc1_lyp_dn;

  *vxc2_upup = xb88_coeff * vx2_b88_upup + clyp_coeff * vc2_lyp_upup;
  *vxc2_dndn = xb88_coeff * vx2_b88_dndn + clyp_coeff * vc2_lyp_dndn;
  *vxc2_updn = xb88_coeff * vx2_b88_updn + clyp_coeff * vc2_lyp_updn;
  *vxc2_dnup = xb88_coeff * vx2_b88_dnup + clyp_coeff * vc2_lyp_dnup;
}
