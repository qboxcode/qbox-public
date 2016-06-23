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
// B3LYPFunctional.C
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
    // not implemented
    assert(false);
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
#if 0 // not implemented
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
      double grx = grx_up + grx_dn;
      double gry = gry_up + gry_dn;
      double grz = grz_up + grz_dn;
      double grad_up = sqrt(grx_up*grx_up + gry_up*gry_up + grz_up*grz_up);
      double grad_dn = sqrt(grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn);
      double grad    = sqrt(grx*grx + gry*gry + grz*grz);
      excpbe_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&exc_up[i],&exc_dn[i],
                &vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],&vxc2_dndn[i],
                &vxc2_updn[i], &vxc2_dnup[i]);
    }
#endif
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
