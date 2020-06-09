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
// VWNFunctional.cpp
//
// LDA Exchange-correlation energy and potential
// S.H.Vosko, L.Wilk, M.Nusair, Can. J. Phys. 58, 1200 (1980)
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include <vector>
#include "VWNFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::setxc(void)
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    // unpolarized
    assert(rho != 0);
    assert(exc != 0);
    assert(vxc1 != 0);
    #pragma omp parallel for
    for ( int ir = 0; ir < _np; ir++ )
    {
      double ex,vx,ec,vc;
      exvwn(rho[ir],ex,vx);
      ecvwn(rho[ir],ec,vc);
      exc[ir] = ex + ec;
      vxc1[ir] = vx + vc;
    }
  }
  else
  {
    // polarized
    assert(rho_up != 0);
    assert(rho_dn != 0);
    assert(exc != 0);
    assert(vxc1_up != 0);
    assert(vxc1_dn != 0);
    #pragma omp parallel for
    for ( int ir = 0; ir < _np; ir++ )
    {
      double ex,vx_up,vx_dn,ec,vc_up,vc_dn;
      exvwn_sp(rho_up[ir],rho_dn[ir],ex,vx_up,vx_dn);
      ecvwn_sp(rho_up[ir],rho_dn[ir],ec,vc_up,vc_dn);
      vxc1_up[ir] = vx_up + vc_up;
      vxc1_dn[ir] = vx_dn + vc_dn;
      exc[ir] = ex + ec;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::exvwn(const double rh, double &ex, double &vx)
{
  x_unpolarized(rh,ex,vx);
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::ecvwn(const double rh, double &ec, double &vc)
{
  c_unpolarized(rh,ec,vc);
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::exvwn_sp(double roe_up, double roe_dn,
  double &ex, double &vx_up, double &vx_dn)
{
  const double fz_prefac = 1.0 / ( cbrt(2.0)*2.0 - 2.0 );
  const double dfz_prefac = (4.0/3.0) * fz_prefac;
  ex = 0.0;
  vx_up = 0.0;
  vx_dn = 0.0;

  if ( roe_up < 0.0 ) roe_up = 0.0;
  if ( roe_dn < 0.0 ) roe_dn = 0.0;
  const double roe = roe_up + roe_dn;

  if ( roe > 0.0 )
  {
    const double zeta = ( roe_up - roe_dn ) / roe;
    const double zp1 = 1.0 + zeta;
    const double zm1 = 1.0 - zeta;
    const double zp1_13 = cbrt(zp1);
    const double zm1_13 = cbrt(zm1);
    const double fz = fz_prefac * ( zp1_13 * zp1 + zm1_13 * zm1 - 2.0 );
    const double dfz = dfz_prefac * ( zp1_13 - zm1_13 );

    double ex_u,vx_u;
    double ex_p,vx_p;
    double a,da;
    x_unpolarized(roe,ex_u,vx_u);
    x_polarized(roe,ex_p,vx_p);
    alpha_c(roe,a,da);

    const double ex_pu = ex_p - ex_u;
    ex = ex_u + fz * ex_pu;
    double vx = vx_u + fz * ( vx_p - vx_u );
    vx_up = vx + ex_pu * ( 1.0 - zeta ) * dfz;
    vx_dn = vx - ex_pu * ( 1.0 + zeta ) * dfz;
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::ecvwn_sp(double roe_up, double roe_dn,
  double &ec, double &vc_up, double &vc_dn)
{
  const double fz_prefac = 1.0 / ( cbrt(2.0)*2.0 - 2.0 );
  const double dfz_prefac = (4.0/3.0) * fz_prefac;
  ec = 0.0;
  vc_up = 0.0;
  vc_dn = 0.0;

  if ( roe_up < 0.0 ) roe_up = 0.0;
  if ( roe_dn < 0.0 ) roe_dn = 0.0;
  const double roe = roe_up + roe_dn;

  if ( roe > 0.0 )
  {
    const double zeta = ( roe_up - roe_dn ) / roe;
    const double zp1 = 1.0 + zeta;
    const double zm1 = 1.0 - zeta;
    const double zp1_13 = cbrt(zp1);
    const double zm1_13 = cbrt(zm1);
    const double fz = fz_prefac * ( zp1_13 * zp1 + zm1_13 * zm1 - 2.0 );
    const double dfz = dfz_prefac * ( zp1_13 - zm1_13 );

    double ec_u,vc_u;
    double ec_p,vc_p;
    double a,da;

    ecvwn(roe,ec_u,vc_u);
    c_polarized(roe,ec_p,vc_p);
    alpha_c(roe,a,da);

    const double zeta3 = zeta*zeta*zeta;
    const double zeta4 = zeta3*zeta;
    a *= (9.0/8.0)*fz_prefac;
    da *= (9.0/8.0)*fz_prefac;
    double ec_pu = ec_p - ec_u - a;

    ec = ec_u + a * fz + ec_pu * fz * zeta4;

    const double vc1 = vc_u + da * fz + ( vc_p - vc_u - da ) * fz * zeta4;
    const double vc2 = a * dfz + ec_pu * zeta3 * ( 4.0*fz + zeta*dfz );

    vc_up = vc1 + ( 1.0 - zeta ) * vc2;
    vc_dn = vc1 - ( 1.0 + zeta ) * vc2;
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::x_unpolarized(const double rh, double &ex, double &vx)
{
  // unpolarized exchange  energy and potential
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711
  const double c3 = -0.610887057711;

  ex = 0.0;
  vx = 0.0;

  if ( rh > 0.0 )
  {
    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    // exchange in Hartree units
    vx = c3 / rs;
    ex = 0.75 * vx;
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::c_unpolarized(const double rh, double &ec, double &vc)
{
  // unpolarized xc energy and potential
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;

  ec = 0.0;
  vc = 0.0;

  if ( rh > 0.0 )
  {
    const double A = 0.0310907;
    const double x0 = -0.10498;
    const double b = 3.72744;
    const double c = 12.9352;
    const double Q = sqrt( 4.0 * c - b * b );
    const double fac1 = 2.0 * b / Q;
    const double fac2 = b * x0 / ( x0 * x0 + b * x0 + c );
    const double fac3 = 2.0 * ( 2.0 * x0 + b ) / Q;

    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    double sqrtrs = sqrt(rs);
    double X = rs + b * sqrtrs + c;
    double fatan = atan( Q / ( 2.0 * sqrtrs + b ) );

    ec = A * ( log( rs / X ) + fac1 * fatan -
               fac2 * ( log( (sqrtrs-x0)*(sqrtrs-x0) / X ) +
                        fac3 * fatan ));

    double t = sqrtrs - x0;
    vc = ec + ( A / 3.0 ) * ( b * sqrtrs * x0 - c * t ) / ( X * t );
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::x_polarized(const double rh, double &ex, double &vx)
{
  // polarized exchange energy and potential
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711
  // c4 = 2**third * c3
  const double c4 = -0.769669463118;

  ex = 0.0;
  vx = 0.0;

  if ( rh > 0.0 )
  {
    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    // Next line : exchange part in Hartree units
    vx = c4 / rs;
    ex = 0.75 * vx;
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::c_polarized(const double rh, double &ec, double &vc)
{
  // polarized correlation energy and potential
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;

  ec = 0.0;
  vc = 0.0;

  if ( rh > 0.0 )
  {
    const double A = 0.01554535;
    const double x0 = -0.32500;
    const double b = 7.06042;
    const double c = 18.0578;
    const double Q = sqrt( 4.0 * c - b * b );
    const double fac1 = 2.0 * b / Q;
    const double fac2 = b * x0 / ( x0 * x0 + b * x0 + c );
    const double fac3 = 2.0 * ( 2.0 * x0 + b ) / Q;

    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    double sqrtrs = sqrt(rs);
    double X = rs + b * sqrtrs + c;
    double fatan = atan( Q / ( 2.0 * sqrtrs + b ) );

    ec = A * ( log( rs / X ) + fac1 * fatan -
               fac2 * ( log( (sqrtrs-x0)*(sqrtrs-x0) / X ) +
                        fac3 * fatan ));

    double t = sqrtrs - x0;
    vc = ec + ( A / 3.0 ) * ( b * sqrtrs * x0 - c * t ) / ( X * t );
  }
}

////////////////////////////////////////////////////////////////////////////////
void VWNFunctional::alpha_c(const double rh, double &a, double &da)
{
  // VWN spin stiffness alpha_c(rh)
  // a = spin stiffness
  // da = d(rh * a)/drh
  // const double third=1.0/3.0;
  // c1 is (3.D0/(4.D0*pi))**third
  const double c1 = 0.6203504908994001;
  // alpha = (4/(9*pi))**third = 0.521061761198
  // const double alpha = 0.521061761198;
  // c2 = -(3/(4*pi)) / alpha = -0.458165293283
  // const double c2 = -0.458165293283;
  // c3 = (4/3) * c2 = -0.610887057711

  a = 0.0;
  da = 0.0;

  if ( rh > 0.0 )
  {
    const double A = 1.0/(6.0*M_PI*M_PI);
    const double x0 = -0.0047584;
    const double b = 1.13107;
    const double c = 13.0045;
    const double Q = sqrt( 4.0 * c - b * b );
    const double fac1 = 2.0 * b / Q;
    const double fac2 = b * x0 / ( x0 * x0 + b * x0 + c );
    const double fac3 = 2.0 * ( 2.0 * x0 + b ) / Q;

    double ro13 = cbrt(rh);
    double rs = c1 / ro13;

    double sqrtrs = sqrt(rs);
    double X = rs + b * sqrtrs + c;
    double fatan = atan( Q / ( 2.0 * sqrtrs + b ) );

    a = A * ( log( rs / X ) + fac1 * fatan -
              fac2 * ( log( (sqrtrs-x0)*(sqrtrs-x0) / X ) +
                       fac3 * fatan ));

    double t = sqrtrs - x0;
    da = a + ( A / 3.0 ) * ( b * sqrtrs * x0 - c * t ) / ( X * t );
  }
}
