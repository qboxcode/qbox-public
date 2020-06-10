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
// SCANFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SCANFunctional.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

SCANFunctional::SCANFunctional(const vector<vector<double> > &rhoe,
  double x_coeff, double c_coeff)
{
  x_coeff_ = x_coeff;
  c_coeff_ = c_coeff;
  _nspin = rhoe.size();
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
  _np = rhoe[0].size();

  if ( _nspin == 1 )
  {
    exc_.resize(_np);
    vxc1_.resize(_np);
    vxc2_.resize(_np);
    vxc3_.resize(_np);
    grad_rho_[0].resize(_np);
    grad_rho_[1].resize(_np);
    grad_rho_[2].resize(_np);
    tau_.resize(_np);
    rho = &rhoe[0][0];
    grad_rho[0] = &grad_rho_[0][0];
    grad_rho[1] = &grad_rho_[1][0];
    grad_rho[2] = &grad_rho_[2][0];
    tau = &tau_[0];
    exc = &exc_[0];
    vxc1 = &vxc1_[0];
    vxc2 = &vxc2_[0];
    vxc3 = &vxc3_[0];
  }
  else
  {
    exc_up_.resize(_np);
    exc_dn_.resize(_np);
    vxc1_up_.resize(_np);
    vxc1_dn_.resize(_np);
    vxc2_upup_.resize(_np);
    vxc2_updn_.resize(_np);
    vxc2_dnup_.resize(_np);
    vxc2_dndn_.resize(_np);
    vxc3_up_.resize(_np);
    vxc3_dn_.resize(_np);
    grad_rho_up_[0].resize(_np);
    grad_rho_up_[1].resize(_np);
    grad_rho_up_[2].resize(_np);
    grad_rho_dn_[0].resize(_np);
    grad_rho_dn_[1].resize(_np);
    grad_rho_dn_[2].resize(_np);
    tau_up_.resize(_np);
    tau_dn_.resize(_np);

    rho_up = &rhoe[0][0];
    rho_dn = &rhoe[1][0];
    grad_rho_up[0] = &grad_rho_up_[0][0];
    grad_rho_up[1] = &grad_rho_up_[1][0];
    grad_rho_up[2] = &grad_rho_up_[2][0];
    grad_rho_dn[0] = &grad_rho_dn_[0][0];
    grad_rho_dn[1] = &grad_rho_dn_[1][0];
    grad_rho_dn[2] = &grad_rho_dn_[2][0];
    tau_up = &tau_up_[0];
    tau_dn = &tau_dn_[0];
    exc_up = &exc_up_[0];
    exc_dn = &exc_dn_[0];
    vxc1_up = &vxc1_up_[0];
    vxc1_dn = &vxc1_dn_[0];
    vxc2_upup = &vxc2_upup_[0];
    vxc2_updn = &vxc2_updn_[0];
    vxc2_dnup = &vxc2_dnup_[0];
    vxc2_dndn = &vxc2_dndn_[0];
    vxc3_up = &vxc3_up_[0];
    vxc3_dn = &vxc3_dn_[0];
  }
}

void SCANFunctional::setxc(void)
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert( rho != 0 );
    assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
    assert( tau != 0 );
    assert( exc != 0 );
    assert( vxc1 != 0 );
    assert( vxc2 != 0 );
    assert( vxc3 != 0 );
    #pragma omp parallel for
    for ( int i = 0; i < _np; i++ )
    {
      double grad = sqrt(grad_rho[0][i]*grad_rho[0][i] +
                         grad_rho[1][i]*grad_rho[1][i] +
                         grad_rho[2][i]*grad_rho[2][i]);
      excSCAN(rho[i], grad, tau[i], &exc[i], &vxc1[i], &vxc2[i], &vxc3[i]);
    }
  }
  else
  {
    assert( rho_up != 0 );
    assert( rho_dn != 0 );
    assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
    assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
    assert( tau_up != 0 );
    assert( tau_dn != 0 );
    assert( exc_up != 0 );
    assert( exc_dn != 0 );
    assert( vxc1_up != 0 );
    assert( vxc1_dn != 0 );
    assert( vxc2_upup != 0 );
    assert( vxc2_updn != 0 );
    assert( vxc2_dnup != 0 );
    assert( vxc2_dndn != 0 );
    assert( vxc3_up != 0 );
    assert( vxc3_dn != 0 );
    #pragma omp parallel for
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
      excSCAN_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,tau_up[i],tau_dn[i],
               &exc_up[i],&exc_dn[i],&vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],
               &vxc2_dndn[i],&vxc2_updn[i],&vxc2_dnup[i],&vxc3_up[i],
               &vxc3_dn[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//  excSCAN: SCAN exchange-correlation
//
//  input:
//    rho:  density
//    grad: abs(grad(rho))
//    tau:  kinetic energy term
//  output:
//    exc: exchange-correlation energy per electron
//    vxc1, vxc2 : quantities such that the total exchange potential is:
//
//      vxc = vxc1 + div ( vxc2 * grad(n) )
//
//    vxc3: additional potential energy term defined as change in energy with
//          respect to kinetic energy tau
//
//  References:
//  [a] S. Jianwei, A. Ruzsinsky, and J.P.Perdew,
//      "Strongly Constrained and Appropriately Normed Semilocal Density
//       Functional", Phys.Rev.Lett. 115, 036402-1, (2015).
//
////////////////////////////////////////////////////////////////////////////////
void SCANFunctional::excSCAN(double rho, double grad, double tau, double *exc,
  double *vxc1, double *vxc2, double *vxc3)
{
  const double pi = M_PI;
  const double hx0 = 1.174;
  const double a1 = 4.9479;
  const double cx1 = 0.667;
  const double cx2 = 0.8;
  const double dx = 1.24;
  const double k1 = 0.065;
  const double mu_AK = 10.0/81.0;
  const double b2 = sqrt(5913.0/405000.0);
  const double b1 = 511.0/(13500.0 * 2.0 * b2);
  const double b3 = 0.5;
  const double b4 = (mu_AK * mu_AK / k1) - 1606.0/18225.0 - b1 * b1;
  const double bc0 = 0.0285764;
  const double bc1 = 0.0889;
  const double bc2 = 0.125541;
  const double cc1 = 0.64;
  const double cc2 = 1.5;
  const double dc = 0.7;
  const double alpha0 = 0.2137;
  const double beta00 = 0.0310907;
  const double beta10 = 7.5957;
  const double beta20 = 3.5876;
  const double beta30 = 1.6382;
  const double beta40 = 0.49294;
  const double gamma = 0.03109069086965489; //gamma = (1-ln2)/pi^2
  const double chi_inf = 0.12802585262625815;
  const double z1 = 0.066724550603149220;
  const double z2 = 0.1;
  const double z3 = 0.1778;

  double tmp0, tmp1, tmp2;

  double rs, rtrs, kF, s, s2, tau_W, tau_unif, XCalpha, oneMalpha;
  double exunif, x, gx, hx1, fx, FXSCAN;
  double dxds, dxdalpha, dgxds, dhx1dx, dfxdalpha;
  double dFXdalpha, dFXds;
  double ex,vx1,vx2,vx3,ec,vc1,vc2,vc3;

  double fc, ec0, ec1, H0, H1, w0, w1;
  double beta1, ginf, A1, t1, g1, g5;
  double ecLDA, decLDAdrs, ecLSDA, decLSDAdrs;

  double dec0drs, dw0drs, dbeta1drs, dw1drs, dA1drs, dt1drs, dH1drs, dec1drs;
  double dt1ds, dg1ds, dec1ds, dginfds, dec0ds;

  double drsdn, dsdn, dalphadn, decdn;
  double dsdgrad, dalphadgrad, decdgrad, dFXdsdsdgrad;
  double dfcdalpha, dalphadtau, decdtau;

  *exc = 0.0;
  *vxc1 = 0.0;
  *vxc2 = 0.0;
  *vxc3 = 0.0;

  if ( rho < 1.e-18  )
  {
    return;
  }

  rs = pow(4.0 * pi * rho / 3.0, -1.0/3.0);
  rtrs = sqrt(rs);
  kF = pow(3.0 * pi * pi * rho, 1.0/3.0);
  s = grad / ( 2.0 * kF * rho );
  s2 = s * s;
  tau_W = grad * grad / (8.0 * rho);
  tau_unif = 0.3 * pow(3.0 * pi * pi, 2.0 / 3.0) * pow(rho, 5.0 / 3.0);
  //!! abs value in next line
  XCalpha = fabs(tau - tau_W) / tau_unif;
  oneMalpha = 1.0 - XCalpha;

  // exchange

  exunif = -3.0 / 4.0 * pow(3.0 * rho / pi, 1.0/3.0);

  // SCAN exchange enhancement factor
  tmp0 = (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
  x = mu_AK * s2 * (1.0 + (b4 * s2 / mu_AK) * exp(-b4 * s2 / mu_AK)) +
      tmp0 * tmp0;
  gx = 1.0 - exp(-a1 / sqrt(s));
  hx1 = 1.0 + k1 - k1/(1.0 + x / k1);
  if ( XCalpha < 1.0 )
  {
    fx = exp(-cx1 * XCalpha / oneMalpha);
  }
  else if ( XCalpha > 1.0 )
  {
    fx = -dx * exp(cx2 / oneMalpha);
  }
  else
  {
    fx = 0.0;
  }

  FXSCAN = gx * (hx1 + fx * (hx0 - hx1));

  //exchange energy
  ex = exunif * FXSCAN;

  // energy done, now the potential
  dxds = 2.0 * s * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK)) + 4.0 * b1 * s *
         (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
  dxdalpha = (2.0 * b3 * oneMalpha * oneMalpha - 1.0) *
             (2.0 * b2 * exp(-b3 * oneMalpha * oneMalpha)) *
             (b1 * s2 + b2 * (oneMalpha) * exp(-b3 * oneMalpha * oneMalpha));
  if ( s == 0 )
  {
    dgxds = 0;
  }
  else
  {
    dgxds = -a1 / (2.0 * pow(s , 1.5)) * exp(-a1 / sqrt(s));
  }
  dhx1dx = (k1 / (k1 + x)) * (k1 / (k1 + x));

  if ( XCalpha < 1.0 )
  {
    dfxdalpha = -cx1 / oneMalpha / oneMalpha * exp(-cx1 * XCalpha / oneMalpha);
  }
  else if ( XCalpha > 1.0 )
  {
    dfxdalpha = -cx2 * dx / oneMalpha / oneMalpha * exp(cx2 / oneMalpha);
  }
  else
  {
    dfxdalpha = 0.0;
  }

  dsdn = -4.0 * s / (3.0 * rho);
  dalphadn = 1.0 / rho * (tau_W / tau_unif - 5.0 / 3.0 * XCalpha);

  dsdgrad = 1.0 / ( 2.0 * kF * rho ) / grad;
  dalphadgrad = -1.0 / (4.0 * rho * tau_unif);
  dalphadtau = 1.0 / tau_unif;

  dFXds = dgxds * (hx1 + fx * (hx0 - hx1)) + gx * (1.0 - fx) * dhx1dx * dxds;
  dFXdalpha = gx * (dhx1dx * dxdalpha + dfxdalpha * (hx0 - hx1) -
                    fx * dhx1dx * dxdalpha);

  if ( s == 0 )
  {
    dFXdsdsdgrad = gx * (1.0 - fx) * dhx1dx *
         (2.0 * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK)) + 4.0 * b1 *
         (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha))) *
         (1.0 / ( 2.0 * kF * rho ) / ( 2.0 * kF * rho ));
  }
  else
  {
    dFXdsdsdgrad = dFXds * dsdgrad;
  }

  //Note: 1/sigma from vx2 distributed to dsdgrad and dalphadgrad to
  //prevent problems with 1/(grad=0)
  vx1 = exunif * ((4.0/3.0) * FXSCAN + rho * (dFXds * dsdn +
                  dFXdalpha * dalphadn));
  vx2 = -rho * exunif * (dFXdsdsdgrad + dFXdalpha * dalphadgrad);
  vx3 = rho * exunif * dFXdalpha * dalphadtau;

  // correlation
  // SCAN Correlation energy
  ecLDA = -bc0/(1.0 + bc1 * rtrs + bc2 * rs);
  ginf = pow(1.0 + 4.0 * chi_inf * s2, -0.25);
  w0 = exp(-ecLDA/bc0) - 1.0;
  H0 = bc0 * log(1.0 + w0 * (1.0 - ginf));
  ec0 = ecLDA + H0;

  beta1 = z1 * ( 1.0 + z2 * rs) / (1.0 + z3 * rs);
  gPW92(alpha0, beta00, beta10, beta20, beta30, beta40, rtrs,
        &ecLSDA, &decLSDAdrs);
  w1 = exp(-ecLSDA/gamma) - 1.0;
  A1 = beta1 / (gamma * w1);
  t1 = pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s / rtrs;
  g1 = pow(1.0 + 4.0 * A1 * t1 * t1, -0.25);
  g5 = g1 * g1 * g1 * g1 * g1;
  H1 = gamma * log(1.0 + w1 * (1.0 - g1));
  ec1 = ecLSDA + H1;
  if ( XCalpha < 1.0 )
  {
    fc = exp(-cc1 * XCalpha / oneMalpha);
  }
  else if ( XCalpha > 1.0 )
  {
    fc = -dc * exp(cc2 / oneMalpha);
  }
  else
  {
    fc = 0.0;
  }

  ec = ec1 + fc * (ec0 - ec1);

  // rs derivatives
  tmp1 = 1.0 + bc1 * rtrs + bc2 * rs;
  decLDAdrs = bc0 * (bc1 / rtrs + 2.0 * bc2) / (2.0 * tmp1 * tmp1);
  dw0drs = -1.0 * (1.0 + w0) / bc0 * decLDAdrs;
  dec0drs = decLDAdrs + (bc0 * (1.0 - ginf))/(1.0 + w0 * (1.0 - ginf)) * dw0drs;

  tmp2 = (1.0 + z3 * rs);
  dbeta1drs = z1 * (z2 - z3) / tmp2 / tmp2;
  dw1drs = - 1.0 * (1.0 + w1) * decLSDAdrs / gamma;
  dA1drs = dbeta1drs / (gamma * w1) - beta1 * dw1drs / (gamma * w1 * w1);
  dt1drs = -1.0 * pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s / (2.0 * rtrs * rtrs *
           rtrs);
  dH1drs = (1.0 - g1) * gamma / (1.0 + w1 * (1.0 - g1)) * dw1drs + w1 * gamma /
           (1.0 + w1 * (1.0 - g1)) * (t1 * t1 * g5 *
           dA1drs + 2.0 * A1 * t1 * g5 * dt1drs);
  dec1drs = decLSDAdrs + dH1drs;

  // s derivatives

  dginfds = -2.0 * chi_inf * s * ginf * ginf * ginf * ginf * ginf;
  dec0ds = -bc0 * w0 * dginfds / (1.0 + w0 * (1.0 - ginf));

  dt1ds = pow(3.0 * pi * pi / 16.0, 1.0/3.0) / rtrs;
  dg1ds = -2.0 * A1 * t1 * g5 * dt1ds;
  dec1ds = -gamma * w1 / (1.0 + w1 * (1.0 - g1)) * dg1ds;

  // alpha derivatives
  if ( XCalpha < 1.0 )
  {
    dfcdalpha = -cc1 / oneMalpha / oneMalpha * exp(-cc1 * XCalpha/ oneMalpha);
  }
  else if ( XCalpha > 1.0 )
  {
    dfcdalpha = -cc2 * dc / oneMalpha / oneMalpha * exp(cc2 / oneMalpha);
  }
  else
  {
    dfcdalpha = 0.0;
  }

  // V1
  drsdn = -rs / (3.0 * rho);
  dsdn = -4.0 * s / (3.0 * rho);
  dalphadn = (tau_W / tau_unif - 5.0 * XCalpha / 3.0) / rho;
  decdn = (dec1drs * drsdn + dec1ds * dsdn) +
          dfcdalpha * dalphadn * (ec0 - ec1) + fc *
          ((dec0drs * drsdn + dec0ds * dsdn) -
          (dec1drs * drsdn + dec1ds * dsdn));
  vc1 = ec + rho * decdn;

  // V2
  dsdgrad = 1.0 / ( 2.0 * kF * rho );
  dalphadgrad = -2.0 * tau_W / (grad * tau_unif);
  decdgrad = dec1ds * dsdgrad + dfcdalpha * dalphadgrad *
    (ec0 - ec1) + fc * ( dec0ds * dsdgrad - dec1ds * dsdgrad );

  // V2 grad=0
  if ( s == 0 )
  {
    double dg1ds2 = -A1 * g5 * dt1ds /
                    (pow(16.0 * rho * rho * rho * rho,1.0/3.0) * rtrs);
    double dec1ds2 = -gamma * w1 / (1.0 + w1 * (1.0 - g1)) * dg1ds2;
    double dginfds2 = -chi_inf * ginf * ginf * ginf * ginf * ginf /
                      pow(3.0 * pi * pi * rho * rho * rho * rho,1.0/3.0);
    double dec0ds2 = -bc0 * w0 * dginfds2 / (1.0 + w0 * (1.0 - ginf));
    double dsdgrad2 = pow(24.0 * pi * pi * rho * rho * rho * rho,-1.0/3.0);
    double dalphadgrad2 = -1.0 / (4.0 * rho * tau_unif);
    double decdgrad2 = dec1ds2 * dsdgrad2 + dfcdalpha * dalphadgrad2 *
    (ec0 - ec1) + fc * (dec0ds2 * dsdgrad2 - dec1ds2 * dsdgrad2);
    vc2 = -rho * decdgrad2;
  }
  else
  {
    vc2 = -rho * decdgrad / grad;
  }

  // V3
  dalphadtau = 1.0 / tau_unif;
  decdtau = dfcdalpha * dalphadtau * (ec0 - ec1);
  vc3 = rho * decdtau;

  // XC
  *exc = x_coeff_ * ex + c_coeff_ * ec;
  *vxc1 = x_coeff_ * vx1 + c_coeff_ * vc1;
  *vxc2 = x_coeff_ * vx2 + c_coeff_ * vc2;
  *vxc3 = x_coeff_ * vx3 + c_coeff_ * vc3;
}

//////////////////////////////////////////////////////////////////////////////
void SCANFunctional::excSCAN_sp(double rho_up, double rho_dn, double grad_up,
    double grad_dn, double grad, double tau_up, double tau_dn, double *exc_up,
    double *exc_dn, double *vxc1_up, double *vxc1_dn, double *vxc2_upup,
    double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup, double *vxc3_up,
    double *vxc3_dn)
{
  const double pi = M_PI;
  const double hx0 = 1.174;
  const double a1 = 4.9479;
  const double cx1 = 0.667;
  const double cx2 = 0.8;
  const double dx = 1.24;
  const double k1 = 0.065;
  const double mu_AK = 10.0/81.0;
  const double b2 = sqrt(5913.0/405000.0);
  const double b1 = 511.0/(13500.0 * 2.0 * b2);
  const double b3 = 0.5;
  const double b4 = (mu_AK * mu_AK / k1) - 1606.0/18225.0 - b1 * b1;
  const double bc0 = 0.0285764;
  const double bc1 = 0.0889;
  const double bc2 = 0.125541;
  const double cc1 = 0.64;
  const double cc2 = 1.5;
  const double dc = 0.7;
  const double GC1 = 2.3631;
  const double alpha0 = 0.2137;
  const double beta00 = 0.0310907;
  const double beta10 = 7.5957;
  const double beta20 = 3.5876;
  const double beta30 = 1.6382;
  const double beta40 = 0.49294;
  const double alpha1 = 0.20548;
  const double beta01 = 0.01554535;
  const double beta11 = 14.1189;
  const double beta21 = 6.1977;
  const double beta31 = 3.3662;
  const double beta41 = 0.62517;
  const double alpha2 = 0.11125;
  const double beta02 = 0.0168869;
  const double beta12 = 10.357;
  const double beta22 = 3.6231;
  const double beta32 = 0.88026;
  const double beta42 = 0.49671;
  const double gamma = 0.03109069086965489; //gamma = (1-ln2)/pi^2
  const double gamma2 = 0.5198420997897463295344212145565;
  const double chi_inf = 0.12802585262625815;
  const double z1 = 0.066724550603149220;
  const double z2 = 0.1;
  const double z3 = 0.1778;

  double tmp0, tmp1, tmp2;

  double rs, rtrs, kF, s, s2, tau_W, tau_unif, ds, XCalpha, oneMalpha;
  double exunif, x, gx, hx1, fx, FXSCAN;
  double dxds, dxdalpha, dgxds, dhx1dx, dfxdalpha;
  double dFXdalpha, dFXds, dFXdsdsdgrad;
  double ex_up,ex_dn,vx1_up,vx1_dn,vx2_up,vx2_dn,vx3_up,vx3_dn,
         ec,vc1_up,vc1_dn,vc2,vc3;
  double rhotot, tautot = tau_up + tau_dn;

  double fc, ec0, ec1, H0, H1, w0, w1, GC;
  double Dx, phi, phi3, zeta, beta1, F, ginf, A1, t1, g1, g5;
  double ecLDA, decLDAdrs, ecLSDA, ecLSDA0, ecLSDA1, ecLSDA2, decLSDAdrs,
         decLSDAdrs0, decLSDAdrs1, decLSDAdrs2;

  double dw0drs, dbeta1drs, dw1drs, dA1drs, dt1drs, dH1drs;
  double dt1ds, dg1ds, dec1ds, dginfds, dH0ds, dec0ds;

  double drsdn, dsdn, dalphadn;
  double dsdgrad, dalphadgrad, decdgrad;
  double dfcdalpha, dalphadtau, decdtau;
  double dalphadzeta, dDxdzeta, dFdzeta, dw1dzeta, decLSDAdzeta,
         dtdzeta, dAdzeta, dgdzeta, dH1dzeta, dphidzeta;
  double dzetadnup, dec1dnup, dH0dnup, dGCdnup, dec0dnup, dfcdnup, decdnup;
  double dzetadndn, dec1dndn, dH0dndn, dGCdndn, dec0dndn, dfcdndn, decdndn;

  *exc_up = 0.0;
  *exc_dn = 0.0;
  *vxc1_up = 0.0;
  *vxc1_dn = 0.0;
  *vxc2_upup = 0.0;
  *vxc2_updn = 0.0;
  *vxc2_dnup = 0.0;
  *vxc2_dndn = 0.0;
  *vxc3_up = 0.0;
  *vxc3_dn = 0.0;

  if ( rho_up < 1.e-18 && rho_dn < 1.e-18  )
  {
    return;
  }

  ex_up = 0.0;
  ex_dn = 0.0;
  vx1_up = 0.0;
  vx1_dn = 0.0;
  vx2_up = 0.0;
  vx2_dn = 0.0;
  vx3_up = 0.0;
  vx3_dn = 0.0;

  //exchange up

  if ( rho_up > 1.e-18 )
  {
    double tworhoup = 2.0 * rho_up;
    double twogradup = 2.0 * grad_up;
    rs = pow(4.0 * pi * tworhoup / 3.0, -1.0/3.0);
    rtrs = sqrt(rs);
    kF = pow(3.0 * pi * pi * tworhoup, 1.0/3.0);
    s = twogradup / ( 2.0 * kF * tworhoup );
    s2 = s * s;
    tau_W = twogradup * twogradup / (8.0 * tworhoup);
    tau_unif = 0.3 * pow(3.0 * pi * pi, 2.0 / 3.0) * pow(tworhoup, 5.0 / 3.0);
    //!! abs value in next line
    XCalpha = fabs(2.0 * tau_up - tau_W) / tau_unif;
    oneMalpha = 1.0 - XCalpha;

    exunif = -3.0 / 4.0 * pow(3.0 * tworhoup / pi, 1.0/3.0);

    // SCAN exchange enhancement factor
    tmp0 = (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
    x = mu_AK * s2 * (1.0 + (b4 * s2 / mu_AK) * exp(-b4 * s2 / mu_AK)) +
        tmp0 * tmp0;
    gx = 1.0 - exp(-a1 / sqrt(s));
    hx1 = 1.0 + k1 - k1/(1.0 + x / k1);
    if ( XCalpha < 1.0 )
    {
      fx = exp(-cx1 * XCalpha / oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      fx = -dx * exp(cx2 / oneMalpha);
    }
    else
    {
      fx = 0.0;
    }

    FXSCAN = gx * (hx1 + fx * (hx0 - hx1));

    //exchange energy
    ex_up = exunif * FXSCAN;

    // energy done, now the potential
    dxds = 2.0 * s * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK)) + 4.0 * b1 * s *
           (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
    dxdalpha = (2.0 * b3 * oneMalpha * oneMalpha - 1.0) *
               (2.0 * b2 * exp(-b3 * oneMalpha * oneMalpha)) *
               (b1 * s2 + b2 * (oneMalpha) * exp(-b3 * oneMalpha * oneMalpha));
    if ( s == 0 )
    {
      dgxds = 0;
    }
    else
    {
      dgxds = -a1 / (2.0 * pow(s , 1.5)) * exp(-a1 / sqrt(s));
    }
    dhx1dx = (k1 / (k1 + x)) * (k1 / (k1 + x));

    if ( XCalpha < 1.0 )
    {
     dfxdalpha = -cx1 / oneMalpha / oneMalpha * exp(-cx1 * XCalpha / oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      dfxdalpha = -cx2 * dx / oneMalpha / oneMalpha * exp(cx2 / oneMalpha);
    }
    else
    {
      dfxdalpha = 0.0;
    }

    dsdn = -4.0 * s / (3.0 * tworhoup);
    dalphadn = 1.0 / tworhoup * (tau_W / tau_unif - 5.0 / 3.0 * XCalpha);

    dsdgrad = 1.0 / ( 2.0 * kF * tworhoup) / twogradup;
    dalphadgrad = -1.0 / (4.0 * tworhoup * tau_unif);
    dalphadtau = 1.0 / tau_unif;

    dFXds = dgxds * (hx1 + fx * (hx0 - hx1)) + gx * (1.0 - fx) * dhx1dx * dxds;
    dFXdalpha = gx * (dhx1dx * dxdalpha + dfxdalpha * (hx0 - hx1) -
                      fx * dhx1dx * dxdalpha);
    if ( s == 0 )
    {
      dFXdsdsdgrad = gx * (1.0 - fx) * dhx1dx *
           (2.0 * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK)) + 4.0 * b1 *
           (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha))) *
           (1.0 / ( 2.0 * kF * tworhoup ) / ( 2.0 * kF * tworhoup ));
    }
    else
    {
      dFXdsdsdgrad = dFXds * dsdgrad;
    }

    //Note: 1/sigma_up from vx2 distributed to dsdgrad and dalphadgrad to
    //prevent problems with 1/(grad=0)
    vx1_up = exunif * ((4.0/3.0) * FXSCAN + tworhoup * (dFXds * dsdn +
                       dFXdalpha * dalphadn));
    vx2_up = -tworhoup * exunif * (dFXdsdsdgrad + dFXdalpha * dalphadgrad);
    vx3_up = tworhoup * exunif * dFXdalpha * dalphadtau;
  }

  // exchange dn

  if ( rho_dn > 1.e-18 )
  {
    double tworhodn = 2.0 * rho_dn;
    double twograddn = 2.0 * grad_dn;
    rs = pow(4.0 * pi * tworhodn / 3.0, -1.0/3.0);
    rtrs = sqrt(rs);
    kF = pow(3.0 * pi * pi * tworhodn, 1.0/3.0);
    s = twograddn / ( 2.0 * kF * tworhodn );
    s2 = s * s;
    tau_W = twograddn * twograddn / (8.0 * tworhodn);
    tau_unif = 0.3 * pow(3.0 * pi * pi, 2.0 / 3.0) * pow(tworhodn, 5.0 / 3.0);
    //!! abs value in next line
    XCalpha = fabs(2.0 * tau_dn - tau_W) / tau_unif;
    oneMalpha = 1.0 - XCalpha;

    exunif = -3.0 / 4.0 * pow(3.0 * tworhodn / pi, 1.0/3.0);

    // SCAN exchange enhancement factor
    tmp0 = (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
    x = mu_AK * s2 * (1.0 + (b4 * s2 / mu_AK) * exp(-b4 * s2 / mu_AK)) +
        tmp0 * tmp0;
    gx = 1.0 - exp(-a1 / sqrt(s));
    hx1 = 1.0 + k1 - k1/(1.0 + x / k1);
    if ( XCalpha < 1.0 )
    {
      fx = exp(-cx1 * XCalpha / oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      fx = -dx * exp(cx2 / oneMalpha);
    }
    else
    {
      fx = 0.0;
    }

    FXSCAN = gx * (hx1 + fx * (hx0 - hx1));

    //exchange energy
    ex_dn = exunif * FXSCAN;

    // energy done, now the potential
    dxds = 2.0 * s * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK)) + 4.0 * b1 * s *
           (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
    dxdalpha = (2.0 * b3 * oneMalpha * oneMalpha - 1.0) *
               (2.0 * b2 * exp(-b3 * oneMalpha * oneMalpha)) *
               (b1 * s2 + b2 * (oneMalpha) * exp(-b3 * oneMalpha * oneMalpha));
    if ( s == 0 )
    {
      dgxds = 0;
    }
    else
    {
      dgxds = -a1 / (2.0 * pow(s , 1.5)) * exp(-a1 / sqrt(s));
    }
    dhx1dx = (k1 / (k1 + x)) * (k1 / (k1 + x));

    if ( XCalpha < 1.0 )
    {
     dfxdalpha = -cx1 / oneMalpha / oneMalpha * exp(-cx1 * XCalpha / oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      dfxdalpha = -cx2 * dx / oneMalpha / oneMalpha * exp(cx2 / oneMalpha);
    }
    else
    {
      dfxdalpha = 0.0;
    }

    dsdn = -4.0 * s / (3.0 * tworhodn);
    dalphadn = 1.0 / tworhodn * (tau_W / tau_unif - 5.0 / 3.0 * XCalpha);

    dsdgrad = 1.0 / ( 2.0 * kF * tworhodn ) / twograddn;
    dalphadgrad = -1.0 / (4.0 * tworhodn * tau_unif);
    dalphadtau = 1.0 / tau_unif;

    dFXds = dgxds * (hx1 + fx * (hx0 - hx1)) + gx * (1.0 - fx) * dhx1dx * dxds;
    dFXdalpha = gx * (dhx1dx * dxdalpha + dfxdalpha * (hx0 - hx1) -
                      fx * dhx1dx * dxdalpha);
    if ( s == 0 )
    {
      dFXdsdsdgrad = gx * (1.0 - fx) * dhx1dx *
           (2.0 * mu_AK * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
           exp(-b4 * s2 / mu_AK)) + 4.0 * b1 *
           (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha))) *
           (1.0 / ( 2.0 * kF * tworhodn ) / ( 2.0 * kF * tworhodn ));
    }
    else
    {
      dFXdsdsdgrad = dFXds * dsdgrad;
    }

    vx1_dn = exunif * ((4.0/3.0) * FXSCAN + tworhodn * (dFXds * dsdn +
                       dFXdalpha * dalphadn));
    vx2_dn = -tworhodn * exunif * (dFXdsdsdgrad + dFXdalpha * dalphadgrad);
    vx3_dn = tworhodn * exunif * dFXdalpha * dalphadtau;
  }

  // set negative densities to positive for correlation part
  if ( rho_up < 0.0 ) rho_up = -rho_up;
  if ( rho_dn < 0.0 ) rho_dn = -rho_dn;
  // set small densities to 0 for correlation part
  if (rho_up < 1.e-18 ) rho_up = 0.0;
  if (rho_dn < 1.e-18 ) rho_dn = 0.0;

  ec = 0.0;
  vc1_up = 0.0;
  vc1_dn = 0.0;
  vc2 = 0.0;
  vc3 = 0.0;

  rhotot = rho_up + rho_dn;

  // correlation
  if ( rhotot > 1.e-18 )
  {
    zeta = (rho_up - rho_dn) / rhotot;
    if ( zeta == 1.0 )
      zeta = 1.0 - 1e-9;
    if ( zeta == -1.0 )
      zeta = -1.0 + 1e-9;
    rs = pow(4.0 * pi * rhotot / 3.0, -1.0/3.0);
    rtrs = sqrt(rs);
    kF = pow(3.0 * pi * pi * rhotot, 1.0/3.0);
    s = grad / ( 2.0 * kF * rhotot );
    s2 = s * s;
    ds = ((pow(1.0 + zeta, 5.0 / 3.0)) + (pow(1.0 - zeta, 5.0 / 3.0))) / 2.0;
    tau_W = grad * grad / (8.0 * rhotot);
    tau_unif = 0.3 * pow(3.0 * pi * pi, 2.0 / 3.0) * pow(rhotot, 5.0 / 3.0) *ds;
    phi = ((pow(1.0 + zeta, 2.0 / 3.0)) + (pow(1.0 - zeta, 2.0 / 3.0))) / 2.0;
    phi3 = phi * phi * phi;
    //!! abs value in next line
    XCalpha = fabs(tautot - tau_W) / tau_unif;
    oneMalpha = 1.0 - XCalpha;

    ecLDA = -bc0/(1.0 + bc1 * rtrs + bc2 * rs);
    ginf = pow(1.0 + 4.0 * chi_inf * s2, -0.25);
    w0 = exp(-ecLDA/bc0) - 1.0;
    H0 = bc0 * log(1.0 + w0 * (1.0 - ginf));
    Dx = ((pow(1.0 + zeta, 4.0 / 3.0)) + (pow(1.0 - zeta, 4.0 / 3.0))) / 2.0;
    GC = (1.0 - GC1 * (Dx - 1.0)) * (1.0 - pow(zeta,12.0));
    ec0 = (ecLDA + H0) * GC;
    beta1 = z1 * ( 1.0 + z2 * rs) / (1.0 + z3 * rs);
    F = (pow(1.0 + zeta,4.0/3.0) + pow(1.0 - zeta,4.0/3.0) - 2.0) / gamma2;
    gPW92(alpha0, beta00, beta10, beta20, beta30, beta40, rtrs,
          &ecLSDA0, &decLSDAdrs0);
    gPW92(alpha1, beta01, beta11, beta21, beta31, beta41, rtrs,
          &ecLSDA1, &decLSDAdrs1);
    gPW92(alpha2, beta02, beta12, beta22, beta32, beta42, rtrs,
          &ecLSDA2, &decLSDAdrs2);
    ecLSDA = ecLSDA0 * (1.0 - F * pow(zeta,4.0)) + ecLSDA1 * F * pow(zeta,4.0)
             - ecLSDA2 * F * (1.0 - pow(zeta,4.0))/(8.0/(9.0 * gamma2));
    w1 = exp(-ecLSDA / (gamma * phi3)) - 1.0;
    A1 = beta1 / (gamma * w1);
    t1 = pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s / (phi * rtrs);
    g1 = pow(1.0 + 4.0 * A1 * t1 * t1, -0.25);
    g5 = g1 * g1 * g1 * g1 * g1;
    H1 = gamma * phi3 * log(1.0 + w1 * (1.0 - g1));
    ec1 = ecLSDA + H1;
    if ( XCalpha < 1.0 )
    {
      fc = exp(-cc1 * XCalpha / oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      fc = -dc * exp(cc2 / oneMalpha);
    }
    else
    {
      fc = 0.0;
    }

    ec = ec1 + fc * (ec0 - ec1);

    // rs derivatives
    tmp1 = 1.0 + bc1 * rtrs + bc2 * rs;
    decLDAdrs = bc0 * (bc1 / rtrs + 2.0 * bc2) / (2.0 * tmp1 * tmp1);
    dw0drs = -1.0 * (1.0 + w0) * decLDAdrs / bc0;

    tmp2 = (1.0 + z3 * rs);
    dbeta1drs = z1 * (z2 - z3) / tmp2 / tmp2;
    decLSDAdrs = decLSDAdrs0 * (1.0 - F * pow(zeta,4.0)) + decLSDAdrs1 * F *
                 pow(zeta,4.0) - decLSDAdrs2 * F * (1.0 - pow(zeta,4.0)) /
                 (8.0 / (9.0 * gamma2));
    dw1drs = - 1.0 * (1.0 + w1) * decLSDAdrs / (gamma * phi3);
    dA1drs = dbeta1drs / (gamma * w1) - beta1 * dw1drs / (gamma * w1 * w1);
    dt1drs = -1.0 * pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s /
             (2.0 * phi * rtrs * rtrs * rtrs);
    dH1drs = (1.0 - g1) * gamma * phi3 * dw1drs / (1.0 + w1 * (1.0 - g1)) +
             w1 * gamma * phi3 / (1.0 + w1 * (1.0 - g1)) *
             (t1 * t1 * g5 * dA1drs + 2.0 * A1 * t1 * g5 * dt1drs);
    //dec1drs = decLSDAdrs + dH1drs;

    // s derivatives

    dginfds = -2.0 * chi_inf * s * ginf * ginf * ginf * ginf * ginf;
    dt1ds = pow(3.0 * pi * pi / 16.0, 1.0/3.0) / phi / rtrs;
    dg1ds = -2.0 * A1 * t1 * g5 * dt1ds;
    dec1ds = -gamma * phi3 * w1 / (1.0 + w1 * (1.0 - g1)) * dg1ds;
    dH0ds = -bc0 * w0 * dginfds / (1.0 + w0 * (1.0 - ginf));
    dec0ds = dH0ds * GC;

    // alpha derivatives
    if ( XCalpha < 1.0 )
    {
      dfcdalpha = -cc1 / oneMalpha / oneMalpha * exp(-cc1 * XCalpha/ oneMalpha);
    }
    else if ( XCalpha > 1.0 )
    {
      dfcdalpha = -cc2 * dc / oneMalpha / oneMalpha * exp(cc2 / oneMalpha);
    }
    else
    {
      dfcdalpha = 0.0;
    }

    // zeta derivatives
    dphidzeta = (pow(1.0 + zeta, -1.0/3.0) - pow(1.0 - zeta, -1.0/3.0)) / 3.0;
    dDxdzeta = 2.0 / 3.0 * (pow(1.0 + zeta, 1.0/3.0) -
               pow(1.0 - zeta, 1.0/3.0));
    dalphadzeta = 5.0 * XCalpha *
                  (pow(1.0 - zeta, 2.0/3.0) - pow(1.0 + zeta, 2.0/3.0)) /
                  (3.0 * (pow(1.0 - zeta, 5.0/3.0) + pow(1.0 + zeta, 5.0/3.0)));
    dFdzeta = 4.0 * (pow(1.0 + zeta, 1.0/3.0) - pow(1.0 - zeta, 1.0/3.0)) /
              (3.0 * gamma2);
    decLSDAdzeta = ecLSDA0 * (-1.0 * dFdzeta * pow(zeta,4.0) -
                   4.0 * F * pow(zeta,3.0)) + ecLSDA1 * (dFdzeta*pow(zeta,4.) +
                   4.0 * F * pow(zeta,3.0)) - ecLSDA2 * (dFdzeta *
                   (1.0 - pow(zeta,4.0)) - 4.0 * F * pow(zeta,3.0)) /
                   (8.0 / (9.0 * gamma2));
    dw1dzeta = (-decLSDAdzeta / (gamma * phi3) + ((3.0 * ecLSDA * dphidzeta) /
                (gamma * phi * phi3))) * exp(-ecLSDA / (gamma * phi3));
    dtdzeta = pow((3.0 * pi * pi) / 16.0, 1.0 / 3.0) * (-s * dphidzeta) /
              (phi * phi * rtrs);
    dAdzeta = (-beta1 * dw1dzeta)/(gamma * w1 * w1);
    dgdzeta = (-dAdzeta * t1 * t1) / pow(1.0 + 4.0 * A1 * t1 * t1, 5.0 / 4.0) +
              (-2.0 * A1 * t1 * dtdzeta) /
              pow(1.0 + 4.0 * A1 * t1 * t1, 5.0 / 4.0);
    dH1dzeta = 3.0 * dphidzeta * H1 / phi + (gamma * phi3 * (dw1dzeta *
               (1.0 - g1))) / (1.0 + w1 * (1.0 - g1)) + (gamma * phi3 *
               (w1 * (-dgdzeta))) / (1.0 + w1 * (1.0 - g1));

    // n derivatives
    drsdn = -rs / (3.0 * rhotot);
    dsdn = -4.0 * s / (3.0 * rhotot);
    dalphadn = (tau_W / tau_unif - 5.0 * XCalpha / 3.0) / rhotot;

    // V1C_up
    dzetadnup = 2.0 * rho_dn / rhotot / rhotot;
    dec1dnup = decLSDAdrs * drsdn + decLSDAdzeta * dzetadnup +
               dH1dzeta * dzetadnup + dH1drs * drsdn + dec1ds * dsdn;
    dH0dnup = bc0 * ((1.0 - ginf) * dw0drs * drsdn - w0 * dginfds * dsdn) /
         (1.0 + w0 * (1.0 - ginf));
    dGCdnup = -dzetadnup * (GC1 * dDxdzeta * (1.0 - pow(zeta,12.0)) +
              12.0 * (1.0 - GC1 * (Dx - 1.0)) * pow(zeta,11.0));
    dec0dnup = (decLDAdrs * drsdn + dH0dnup) * GC + (ecLDA + H0) * dGCdnup;
    dfcdnup = dfcdalpha * (dalphadn + dalphadzeta * dzetadnup);
    decdnup = dec1dnup + dfcdnup * (ec0 - ec1) + fc * (dec0dnup - dec1dnup);
    vc1_up = ec + rhotot * decdnup;

    // V1C_dn
    dzetadndn = -2.0 * rho_up / rhotot / rhotot;
    dec1dndn = decLSDAdrs * drsdn + decLSDAdzeta * dzetadndn +
               dH1dzeta * dzetadndn + dH1drs * drsdn + dec1ds * dsdn;
    dH0dndn = bc0 * ((1.0 - ginf) * dw0drs * drsdn - w0 * dginfds * dsdn) /
         (1.0 + w0 * (1.0 - ginf));
    dGCdndn = -dzetadndn * (GC1 * dDxdzeta * (1.0 - pow(zeta,12.0)) +
              12.0 * (1.0 - GC1 * (Dx - 1.0)) * pow(zeta,11.0));
    dec0dndn = (decLDAdrs * drsdn + dH0dndn) * GC + (ecLDA + H0) * dGCdndn;
    dfcdndn = dfcdalpha * (dalphadn + dalphadzeta * dzetadndn);
    decdndn = dec1dndn + dfcdndn * (ec0 - ec1) + fc * (dec0dndn - dec1dndn);
    vc1_dn = ec + rhotot * decdndn;

    // VC2
    dsdgrad = s / grad;
    dalphadgrad = -2.0 * tau_W / (grad * tau_unif);
    decdgrad = dec1ds * dsdgrad + dfcdalpha * dalphadgrad *
               (ec0 - ec1) + fc * ( dec0ds * dsdgrad - dec1ds * dsdgrad );

    // V2 grad=0
    if ( s == 0 )
    {
      double dg1ds2 = -A1 * g5 * dt1ds /
                      (pow(16.0 * pow(rhotot,4.0),1.0/3.0) * phi * rtrs);
      double dec1ds2 = -gamma * w1 / (1.0 + w1 * (1.0 - g1)) * dg1ds2;
      double dginfds2 = -chi_inf * ginf * ginf * ginf * ginf * ginf /
                        pow(3.0 * pi * pi * pow(rhotot,4.0),1.0/3.0);
      double dec0ds2 = -bc0 * w0 * dginfds2 * GC / (1.0 + w0 * (1.0 - ginf));
      double dsdgrad2 = pow(24.0 * pi * pi * pow(rhotot,5.0),-1.0/3.0);
      double dalphadgrad2 = -1.0 / (4.0 * rhotot * tau_unif);
      double decdgrad2 = dec1ds2 * dsdgrad2 + dfcdalpha * dalphadgrad2 *
      (ec0 - ec1) + fc * (dec0ds2 * dsdgrad2 - dec1ds2 * dsdgrad2);
      vc2 = -rhotot * decdgrad2;
    }
    else
    {
      vc2 = -rhotot * decdgrad / grad;
    }

    // VC3
    dalphadtau = 1.0 / tau_unif;
    decdtau = dfcdalpha * dalphadtau * (ec0 - ec1);
    vc3 = rhotot * decdtau;
  }
  *exc_up = x_coeff_ * ex_up + c_coeff_ * ec;
  *exc_dn = x_coeff_ * ex_dn + c_coeff_ * ec;
  *vxc1_up = x_coeff_ * vx1_up + c_coeff_ * vc1_up;
  *vxc1_dn = x_coeff_ * vx1_dn + c_coeff_ * vc1_dn;
  *vxc2_upup = x_coeff_ * 2.0 * vx2_up + c_coeff_ * vc2;
  *vxc2_dndn = x_coeff_ * 2.0 * vx2_dn + c_coeff_ * vc2;
  *vxc2_updn = c_coeff_ * vc2;
  *vxc2_dnup = c_coeff_ * vc2;
  *vxc3_up = x_coeff_ * vx3_up + c_coeff_ * vc3;
  *vxc3_dn = x_coeff_ * vx3_dn + c_coeff_ * vc3;
}

////////////////////////////////////////////////////////////////////////////////
//
//  gPW92.c: Interpolate LSD correlation energy
//  as given by (10) of Perdew & Wang, Phys Rev B45 13244 (1992)
//  Translated into C Dec 9, 1996
//
////////////////////////////////////////////////////////////////////////////////

void SCANFunctional::gPW92(double alpha, double beta0, double beta1,
  double beta2, double beta3, double beta4, double rtrs, double *gg,
  double *dgdrs)
{
  double q0,q1,q2,q3;
  q0 = -2.0 * beta0 * ( 1.0 + alpha * rtrs * rtrs );
  q1 = 2.0 * beta0 * rtrs * ( beta1 + rtrs * ( beta2 + rtrs *
       ( beta3 + rtrs * beta4 ) ) );
  q2 = log ( 1.0 + 1.0 / q1 );
  q3 = beta0 * ( beta1 / rtrs + 2.0 * beta2 + rtrs *
       ( 3.0 * beta3 + 4.0 * beta4 * rtrs ));
  *gg = q0 * q2;
  *dgdrs = -2.0 * alpha * beta0 * q2 - q0 * q3 / ( q1 * ( 1.0 + q1 ));
}
