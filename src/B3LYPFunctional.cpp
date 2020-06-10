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
// B3LYPFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include "B3LYPFunctional.h"
#include "BLYPFunctional.h"
#include "VWNFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
B3LYPFunctional::B3LYPFunctional(const vector<vector<double> > &rhoe)
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
void B3LYPFunctional::setxc(void)
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
      excb3lyp(rho[i],grad,&exc[i],&vxc1[i],&vxc2[i]);
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

      excb3lyp_sp(rho_up[i],rho_dn[i],grad_up2,grad_dn2,grad_up_grad_dn,
                  &exc_up[i],&exc_dn[i],&vxc1_up[i],&vxc1_dn[i],
                  &vxc2_upup[i],&vxc2_dndn[i],
                  &vxc2_updn[i],&vxc2_dnup[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void B3LYPFunctional::excb3lyp(double rho, double grad,
  double *exc, double *vxc1, double *vxc2)
{
  // B3LYP unpolarized
  // Coefficients of the B3LYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // 0.2 ExHF + 0.80 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  const double xlda_coeff = 0.80; // Slater exchange
  const double clda_coeff = 0.19; // LSDA correlation
  const double xb88_coeff = 0.72; // Becke88 exchange gradient correction
  const double clyp_coeff = 1.0 - clda_coeff;

  double ex_lda,vx_lda,ec_lda,vc_lda;
  double ex_b88,vx1_b88,vx2_b88;
  double ec_lyp,vc1_lyp,vc2_lyp;

  VWNFunctional::exvwn(rho,ex_lda,vx_lda);
  VWNFunctional::ecvwn(rho,ec_lda,vc_lda);
  BLYPFunctional::exb88(rho,grad,&ex_b88,&vx1_b88,&vx2_b88);
  BLYPFunctional::eclyp(rho,grad,&ec_lyp,&vc1_lyp,&vc2_lyp);

  const double dex_b88 = ex_b88 - ex_lda;
  const double dvx1_b88 = vx1_b88 - vx_lda;
  const double dvx2_b88 = vx2_b88;

  *exc = xlda_coeff  * ex_lda +
         clda_coeff  * ec_lda +
         xb88_coeff  * dex_b88 +
         clyp_coeff  * ec_lyp;

  *vxc1 = xlda_coeff * vx_lda +
          clda_coeff * vc_lda +
          xb88_coeff * dvx1_b88 +
          clyp_coeff * vc1_lyp;

  *vxc2 = xb88_coeff * dvx2_b88 +
          clyp_coeff * vc2_lyp;
}

////////////////////////////////////////////////////////////////////////////////
void B3LYPFunctional::excb3lyp_sp(double rho_up, double rho_dn,
  double grad_up2, double grad_dn2, double grad_up_grad_dn,
  double *exc_up, double *exc_dn, double *vxc1_up, double *vxc1_dn,
  double *vxc2_upup, double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup)
{
  // B3LYP spin-polarized
  // Coefficients of the B3LYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // 0.2 ExHF + 0.80 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  const double xlda_coeff = 0.80; // Slater exchange
  const double clda_coeff = 0.19; // LSDA correlation
  const double xb88_coeff = 0.72; // Becke88 exchange gradient correction
  const double clyp_coeff = 1.0 - clda_coeff;

  double ex_lda,vx_lda_up,vx_lda_dn;
  double ec_lda,vc_lda_up,vc_lda_dn;
  double ex_b88_up,ex_b88_dn,vx1_b88_up,vx1_b88_dn,
         vx2_b88_upup,vx2_b88_dndn,vx2_b88_updn,vx2_b88_dnup;
  double ec_lyp_up,ec_lyp_dn,vc1_lyp_up,vc1_lyp_dn,
         vc2_lyp_upup,vc2_lyp_dndn,vc2_lyp_updn,vc2_lyp_dnup;

  VWNFunctional::exvwn_sp(rho_up,rho_dn,ex_lda,vx_lda_up,vx_lda_dn);
  VWNFunctional::ecvwn_sp(rho_up,rho_dn,ec_lda,vc_lda_up,vc_lda_dn);
  BLYPFunctional::exb88_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ex_b88_up,&ex_b88_dn,&vx1_b88_up,&vx1_b88_dn,
    &vx2_b88_upup,&vx2_b88_dndn,&vx2_b88_updn,&vx2_b88_dnup);
  BLYPFunctional::eclyp_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ec_lyp_up,&ec_lyp_dn,&vc1_lyp_up,&vc1_lyp_dn,
    &vc2_lyp_upup,&vc2_lyp_dndn,&vc2_lyp_updn,&vc2_lyp_dnup);

  const double dex_b88_up = ex_b88_up - ex_lda;
  const double dex_b88_dn = ex_b88_dn - ex_lda;
  const double dvx1_b88_up = vx1_b88_up - vx_lda_up;
  const double dvx1_b88_dn = vx1_b88_dn - vx_lda_dn;
  const double dvx2_b88_upup = vx2_b88_upup;
  const double dvx2_b88_dndn = vx2_b88_dndn;
  const double dvx2_b88_updn = vx2_b88_updn;
  const double dvx2_b88_dnup = vx2_b88_dnup;

  *exc_up = xlda_coeff  * ex_lda +
            clda_coeff  * ec_lda +
            xb88_coeff  * dex_b88_up +
            clyp_coeff  * ec_lyp_up;

  *exc_dn = xlda_coeff  * ex_lda +
            clda_coeff  * ec_lda +
            xb88_coeff  * dex_b88_dn +
            clyp_coeff  * ec_lyp_dn;

  *vxc1_up = xlda_coeff * vx_lda_up +
             clda_coeff * vc_lda_up +
             xb88_coeff * dvx1_b88_up +
             clyp_coeff * vc1_lyp_up;

  *vxc1_dn = xlda_coeff * vx_lda_dn +
             clda_coeff * vc_lda_dn +
             xb88_coeff * dvx1_b88_dn +
             clyp_coeff * vc1_lyp_dn;

  *vxc2_upup = xb88_coeff * dvx2_b88_upup + clyp_coeff * vc2_lyp_upup;
  *vxc2_dndn = xb88_coeff * dvx2_b88_dndn + clyp_coeff * vc2_lyp_dndn;
  *vxc2_updn = xb88_coeff * dvx2_b88_updn + clyp_coeff * vc2_lyp_updn;
  *vxc2_dnup = xb88_coeff * dvx2_b88_dnup + clyp_coeff * vc2_lyp_dnup;
}
