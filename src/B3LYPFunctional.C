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

  *exc = 0.0;
  *vxc1 = 0.0;
  *vxc2 = 0.0;

  if ( rho < 1.e-18  )
  {
    return;
  }

  // LDA correlation
  // Perdew-Zunger parametrization of Ceperley-Alder data
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711
  const double c3 = -0.610887057711;

  const double A  =  0.0311;
  const double B  = -0.048;
  const double b1 =  1.0529;
  const double b2 =  0.3334;
  const double G  = -0.1423;

  // C from the PZ paper: const double C  =  0.0020;
  // D from the PZ paper: const double D  = -0.0116;
  // C and D by matching Ec and Vc at rs=1
  const double D = G / ( 1.0 + b1 + b2 ) - B;
  const double C = -A - D - G * ( (b1/2.0 + b2) / ((1.0+b1+b2)*(1.0+b1+b2)));

  double ro13 = cbrt(rho);
  double rs = c1 / ro13;

  double ec_lda=0.0,vc_lda=0.0;

  // Next line : exchange in Hartree units
  double vx_lda = c3 / rs;
  double ex_lda = 0.75 * vx_lda;

  // Next lines : Ceperley & Alder correlation (Zunger & Perdew)
  if ( rs < 1.0 )
  {
    double logrs = log(rs);
    ec_lda = A * logrs + B + C * rs * logrs + D * rs;
    vc_lda = A * logrs + ( B - A / 3.0 ) +
            (2.0/3.0) * C * rs * logrs +
            ( ( 2.0 * D - C ) / 3.0 ) * rs;
  }
  else
  {
    double sqrtrs = sqrt(rs);
    double den = 1.0 + b1 * sqrtrs + b2 * rs;
    ec_lda = G / den;
    vc_lda = ec_lda * ( 1.0 + (7.0/6.0) * b1 * sqrtrs +
                      (4.0/3.0) * b2 * rs ) / den;
  }

  // Becke88 exchange: A.D.Becke, Phys.Rev. B38, 3098 (1988)
  // Becke88 exchange constants
  const double beta=0.0042;
  //const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double axa = -0.9305257363490999; /* -1.5*pow(3.0/(4*pi),third) */

  const double rha = 0.5 * rho;
  const double grada = 0.5 * grad;

  const double rha13 = pow ( rha, 1.0/3.0 );
  const double rha43 = rha * rha13;
  const double xa = grada / rha43;
  const double xa2 = xa*xa;
  const double asinhxa = asinh(xa);
  const double frac = 1.0 / ( 1.0 + 6.0 * beta * xa * asinhxa );
  const double ga = axa - beta * xa2 * frac;

  const double ex_b88 = rha13 * ga;

  // Becke88 GGA exchange correction
  const double dex_b88 = ex_b88 - ex_lda;

      // Becke88 potential
  const double gpa = ( 6.0*beta*beta*xa2 * ( xa/sqrt(xa2+1.0) - asinhxa )
                     - 2.0*beta*xa ) * frac*frac;
  const double vx1_b88 = rha13 * (4.0/3.0) * ( ga - xa * gpa );
  const double vx2_b88 = - 0.5 * gpa / grada;

  // Becke88 GGA exchange correction potential
  const double dvx1_b88 = vx1_b88 - vx_lda;
  const double dvx2_b88 = vx2_b88;

  //------------------------------------------------------------
  // LYP correlation
  // Phys. Rev. B 37, 785 (1988).
  // LYP constants
  const double a = 0.04918;
  const double b = 0.132;
  const double ab36 = a * b / 36.0;
  const double c = 0.2533;
  const double c_third = c / 3.0;
  const double d = 0.349;
  const double d_third = d / 3.0;
  const double cf = 2.87123400018819; /* (3/10)*pow(3*pi*pi,2/3) */
  const double cfb = cf * b;

  // next lines specialized to the unpolarized case
  const double rh13 = pow ( rho, 1.0/3.0 );
  const double rhm13 = 1.0 / rh13;
  const double rhm43 = rhm13 / rho;
  const double e = exp ( - c * rhm13 );
  const double num = 1.0 + cfb * e;
  const double den = 1.0 + d * rhm13;
  const double deninv = 1.0 / den;
  const double cfrac = num * deninv;

  const double delta = rhm13 * ( c + d * deninv );
  const double rhm53 = rhm43 * rhm13;
  const double t1 = e * deninv;
  const double t2 = rhm53;
  const double t3 = 6.0 + 14.0 * delta;

  const double g = ab36 * t1 * t2 * t3;

  /* next line, ec is the energy density, hence divide the energy by rho */
  const double ec_lyp = - a * cfrac + 0.25 * g * grad * grad / rho;

  /* energy done, now the potential */
  const double de = c_third * rhm43 * e;
  const double dnum = cfb * de;
  const double dden = - d_third * rhm43;
  const double dfrac = ( dnum * den - dden * num ) * deninv * deninv;

  const double ddelta = - (1.0/3.0) * rhm43 * ( c + d * deninv ) -
           rhm13 * d * dden * deninv * deninv;
  const double dt1 = de * deninv - e * dden * deninv * deninv;
  const double dt2 = - (5.0/3.0) * rhm53/rho;
  const double dt3 = 14.0 * ddelta;

  const double dg = ab36 * ( dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3 );

  const double vc1_lyp = - a * (cfrac + rho * dfrac) + 0.25 * dg * grad * grad;
  const double vc2_lyp = -0.5 * g;

  // Coefficients of the B3LYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Perdew-Zunger parametrization of Ceperley-Alder data
  // 0.2 ExHF + 0.80 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  const double xlda_coeff = 0.80; // Slater exchange
  const double clda_coeff = 0.19; // LSDA correlation
  const double xb88_coeff = 0.72; // Becke88 exchange gradient correction
  const double clyp_coeff = 1.0 - clda_coeff;

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
