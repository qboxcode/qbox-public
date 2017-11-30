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
// SCANFunctional.C
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
    grad_rho_up_[0].resize(_np);
    grad_rho_up_[1].resize(_np);
    grad_rho_up_[2].resize(_np);
    grad_rho_dn_[0].resize(_np);
    grad_rho_dn_[1].resize(_np);
    grad_rho_dn_[2].resize(_np);

    rho_up = &rhoe[0][0];
    rho_dn = &rhoe[1][0];
    grad_rho_up[0] = &grad_rho_up_[0][0];
    grad_rho_up[1] = &grad_rho_up_[1][0];
    grad_rho_up[2] = &grad_rho_up_[2][0];
    grad_rho_dn[0] = &grad_rho_dn_[0][0];
    grad_rho_dn[1] = &grad_rho_dn_[1][0];
    grad_rho_dn[2] = &grad_rho_dn_[2][0];
    exc_up = &exc_up_[0];
    exc_dn = &exc_dn_[0];
    vxc1_up = &vxc1_up_[0];
    vxc1_dn = &vxc1_dn_[0];
    vxc2_upup = &vxc2_upup_[0];
    vxc2_updn = &vxc2_updn_[0];
    vxc2_dnup = &vxc2_dnup_[0];
    vxc2_dndn = &vxc2_dndn_[0];
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
                         grad_rho[2][i]*grad_rho[2][i] );
      excSCAN(rho[i], grad, tau[i], &exc[i], &vxc1[i], &vxc2[i], &vxc3[i]);
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
      excSCAN_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&exc_up[i],&exc_dn[i],
                &vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],&vxc2_dndn[i],
                &vxc2_updn[i], &vxc2_dnup[i]);
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
  const double beta00 = 0.031097;
  const double beta10 = 7.5957;
  const double beta20 = 3.5876;
  const double beta30 = 1.6382;
  const double beta40 = 0.49294;
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double chi_inf = 0.128026;

  double tmp0, tmp1, tmp2;

  double rs, rtrs, kF, s, s2, tau_W, tau_unif, XCalpha, oneMalpha;
  double exunif, x, gx, hx1, fx, fxSCAN;
  double dxds, dxdalpha, dgxds, dhx1ds, dhx1dalpha, dfxdalpha;
  double fXalpha, fXs;
  double ex,vx1,vx2,vx3,ec,vc1,vc2,vc3;

  double fc, ec0, ec1, H0, H1, w0, w1;
  double beta1, ginf, A1, t1, g1;
  double ecLDA, decLDAdrs, ecLSDA, decLSDAdrs;

  double dec0drs, dw0drs, dbeta1drs, dw1drs, dA1drs, dt1drs, dH1drs, dec1drs;
  double dt1ds, dg1ds, dec1ds, dginfds, dec0ds;

  double drsdn, dsdn, dalphadn, decdn;
  double dsdgrad, dalphadgrad, decdgrad;
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
  s  = grad / ( 2.0 * kF * rho );
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
  else if (XCalpha > 1.0)
  {
    fx = -dx * exp(cx2 / oneMalpha);
  }
  else
  {
    fx = 0.0;
  }
  fxSCAN = gx * (hx1 + fx * (hx0 - hx1));

  //exchange energy
  ex = exunif * fxSCAN;

  // energy done, now the potential
  dxds = 2.0 * mu_AK * s * (1.0 + 2.0 * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK) - (b4 * s2/ mu_AK) * (b4 * s2/ mu_AK) *
         exp(-b4 * s2 / mu_AK)) + 4.0 * b1 * s *
         (b1 * s2 + b2 * oneMalpha * exp(-b3 * oneMalpha * oneMalpha));
  dxdalpha = (2.0 * b3 * oneMalpha * oneMalpha - 1.0) *
             (2.0 * b2 * exp(-b3 * oneMalpha * oneMalpha)) *
             (b1 * s2 + b2 * (oneMalpha) * exp(-b3 * oneMalpha * oneMalpha));
  dgxds= -a1/(2.0 * pow(s , 1.5)) * exp(-a1 / sqrt(s));
  dhx1ds = k1 * k1 * dxds / (k1 + x) / (k1 + x);
  dhx1dalpha = k1 * k1 * dxdalpha / (k1 + x) / (k1 + x);

  if (XCalpha < 1.0)
  {
    dfxdalpha = -cx1 / oneMalpha / oneMalpha * exp(-cx1 * XCalpha / oneMalpha);
  }
  else if (XCalpha > 1.0)
  {
    dfxdalpha = -cx2 * dx / oneMalpha / oneMalpha * exp(cx2 / oneMalpha);
  }
  else
  {
    dfxdalpha = 0.0;
  }

  fXs = dgxds * (hx1 + fx * (hx0 - hx1)) + gx * dhx1ds * (1.0 - fx);
  fXalpha = gx * (dhx1dalpha + dfxdalpha * (hx0 - hx1) - fx * dhx1dalpha);

  vx1 = (4.0/3.0) * exunif * (fxSCAN - s * fXs +
        fXalpha * (3.0 * tau_W/(4.0 * tau_unif) - 5.0 * XCalpha / 4.0));
  vx2 = -rho / (grad*grad) * exunif *
        (s * fXs - 2.0 * fXalpha * tau_W/tau_unif);
  vx3 = exunif * (rho / tau_unif) * fXalpha;


  // correlation
  // SCAN Correlation energy
  ecLDA = -bc0/(1.0 + bc1 * rtrs + bc2 * rs);
  ginf = pow(1.0 + 4.0 * chi_inf * s2, -0.25);
  w0 = exp(-ecLDA/bc1) - 1.0;
  H0 = bc1 * log(1.0 + w0 * (1.0 - ginf));
  ec0 = ecLDA + H0;

  beta1 = 0.066725 * ( 1.0 + 0.1 * rs) / (1.0 + 0.1778 * rs);
  gPW92(alpha0, beta00, beta10, beta20, beta30, beta40, rtrs,
        &ecLSDA, &decLSDAdrs);
  w1 = exp(-ecLSDA/gamma) - 1.0;
  A1 = beta1 / (gamma * w1);
  t1 = pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s / rtrs;
  g1 = pow(1.0 + 4.0 * A1 * t1 * t1, -0.25);
  H1 = gamma * log(1.0 + w1 * (1.0 - g1));
  ec1 = ecLSDA + H1;
  if ( XCalpha < 1.0)
  {
    fc = exp(-cc1 * XCalpha / oneMalpha);
  }
  else if (XCalpha > 1.0)
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
  dw0drs = -1.0 * (1.0 + w0) / bc1 * decLDAdrs;
  dec0drs = decLDAdrs + (bc1 * (1.0 - ginf))/(1.0 + w0 * (1.0 - ginf)) * dw0drs;

  tmp2 = (1.0 + 0.1778 * rs);
  dbeta1drs = 0.066725 * (0.1 - 0.1778) / tmp2 / tmp2;
  dw1drs = - 1.0 * (1.0 + w1) * decLSDAdrs / gamma;
  dA1drs = dbeta1drs / (gamma * w1) - beta1 * dw1drs / (gamma * w1 * w1);
  dt1drs = -1.0 * pow(3.0 * pi * pi / 16.0, 1.0/3.0) * s / (2.0 * rtrs * rtrs *
           rtrs);
  dH1drs = (1.0 - g1) * gamma / (1.0 + w1 * (1.0 - g1)) * dw1drs + w1 * gamma /
           (1.0 + w1 * (1.0 - g1)) * (t1 * t1 * g1 * g1 * g1 * g1 * g1 *
           dA1drs + 2.0 * A1 * t1 * g1 * g1 * g1 * g1 * g1 * dt1drs);
  dec1drs = decLSDAdrs + dH1drs;

  // s derivatives

  dginfds = -2.0 * chi_inf * s * ginf * ginf * ginf * ginf * ginf;
  dec0ds = -bc1 * w0 * dginfds / (1.0 + w0 * (1.0 - ginf));

  dt1ds = pow(3.0 * pi * pi / 16.0, 1.0/3.0) / rtrs;
  dg1ds = -2.0 * A1 * t1 * g1 * g1 * g1 * g1 * g1 * dt1ds;
  dec1ds = -gamma * w1 / (1.0 + w1 * (1.0 - g1)) * dg1ds;

  // alpha derivatives
  if ( XCalpha < 1.0)
  {
    dfcdalpha = -cc1 / oneMalpha / oneMalpha * exp(-cc1 * XCalpha/ oneMalpha);
  }
  else if (XCalpha > 1.0)
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
  dsdgrad = s / grad;
  dalphadgrad = -2.0 * tau_W / (grad * tau_unif);
  decdgrad = dec1ds * dsdgrad + dfcdalpha * dalphadgrad *
    (ec0 - ec1) + fc * ( dec0ds * dsdgrad - dec1ds * dsdgrad );
  vc2 = -rho * decdgrad / grad;

  // V3
  dalphadtau = 1.0 / tau_unif;
  decdtau = dfcdalpha * dalphadtau * (ec0 - ec1);
  vc3 = - rho * decdtau;

  // XC
  *exc = x_coeff_ * ex + c_coeff_ * ec;
  *vxc1 = x_coeff_ * vx1 + c_coeff_ * vc1;
  *vxc2 = x_coeff_ * vx2 + c_coeff_ * vc2;
  *vxc3 = x_coeff_ * vx3 + c_coeff_ * vc3;
}

//////////////////////////////////////////////////////////////////////////////

void SCANFunctional::excSCAN_sp(double rho_up, double rho_dn,
  double grad_up, double grad_dn, double grad, double *exc_up, double *exc_dn,
  double *vxc1_up, double *vxc1_dn, double *vxc2_upup, double *vxc2_dndn,
  double *vxc2_updn, double *vxc2_dnup)
{
  const double third  = 1.0 / 3.0;
  const double third2 =  2.0 / 3.0;
  const double third4 =  4.0 / 3.0;
  const double sixthm = -1.0 / 6.0;
  const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double um = 0.2195149727645171;
  const double uk = 0.804;
  const double ul = um / uk;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
  const double alpha = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
  const double seven_sixth  =  7.0 / 6.0;
  const double four_over_pi = 1.27323954473516;
  const double gam = 0.5198420997897463; /* gam = 2^(4/3) - 2 */
  const double fzz = 8.0 / ( 9.0 * gam );
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double bet = 0.06672455060314922; /* see [a] (4) */
  const double delt = bet / gamma;
  const double eta = 1.e-12; // small number to avoid blowup as |zeta|->1

  double eu,eurs,ep,eprs,alfm,alfrsm;
  double ex_up,ex_dn,vx1_up,vx1_dn,vx2_up,vx2_dn,ec,vc1_up,vc1_dn,vc2;

  *exc_up = 0.0;
  *exc_dn = 0.0;
  *vxc1_up = 0.0;
  *vxc1_dn = 0.0;
  *vxc2_upup = 0.0;
  *vxc2_updn = 0.0;
  *vxc2_dnup = 0.0;
  *vxc2_dndn = 0.0;

  if ( rho_up < 1.e-18 && rho_dn < 1.e-18  )
  {
    return;
  }

  /* exchange up */

  ex_up = 0.0;
  vx1_up = 0.0;
  vx2_up = 0.0;
  if ( rho_up > 1.e-18 )
  {
    double tworho = 2.0 * rho_up;
    double gr = 2.0 * grad_up;

    double rh13 = pow ( tworho, third );
    /* LDA exchange energy density */
    double exunif = ax * rh13;
    /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
    double fk = pi32third * rh13;
    double s  = gr / ( 2.0 * fk * tworho );
    /* SCAN enhancement factor */
    double s2 = s * s;
    double p0 = 1.0 + ul * s2;
    double fxSCAN = 1.0 + uk - uk / p0;
    ex_up = exunif * fxSCAN;
    /* energy done, now the potential */
    /* find first derivative of Fx w.r.t the variable s. */
    /* fs = (1/s) * d Fx / d s */
    double fs = 2.0 * uk * ul / ( p0 * p0 );
    vx1_up = third4 * exunif * ( fxSCAN - s2 * fs );
    vx2_up = - exunif * fs / ( tworho * 4.0 * fk * fk );
  }

  /* exchange dn */

  ex_dn = 0.0;
  vx1_dn = 0.0;
  vx2_dn = 0.0;
  if ( rho_dn > 1.e-18 )
  {
    double tworho = 2.0 * rho_dn;
    double gr = 2.0 * grad_dn;

    double rh13 = pow ( tworho, third );
    /* LDA exchange energy density */
    double exunif = ax * rh13;
    /* Fermi wavevector  kF = ( 3 * pi^2 n )^(1/3) */
    double fk = pi32third * rh13;
    double s  = gr / ( 2.0 * fk * tworho );
    /* SCAN enhancement factor */
    double s2 = s * s;
    double p0 = 1.0 + ul * s2;
    double fxSCAN = 1.0 + uk - uk / p0;
    ex_dn = exunif * fxSCAN;
    /* energy done, now the potential */
    /* find first derivative of Fx w.r.t the variable s. */
    /* fs = (1/s) * d Fx / d s */
    double fs = 2.0 * uk * ul / ( p0 * p0 );
    vx1_dn = third4 * exunif * ( fxSCAN - s2 * fs );
    vx2_dn = - exunif * fs / ( tworho * 4.0 * fk * fk );
  }

  /* set negative densities to 0 for correlation part */

  if ( rho_up < 1.e-18 ) rho_up=0.0;
  if ( rho_dn < 1.e-18 ) rho_dn=0.0;

  /* correlation */

  // Find LSD contributions, using [c] (10) and Table I of [c].
  // eu = unpolarized LSD correlation energy
  // eurs = d eu / d rs
  // ep = fully polarized LSD correlation energy
  // eprs = d ep / d rs
  // alfm = - spin stiffness, [c] (3)
  // alfrsm = -d alpha / d rs
  // f = spin-scaling factor from [c] (9)
  // construct ec, using [c] (8)

  double rhotot = rho_up + rho_dn;

  double rh13 = pow ( rhotot, third );
  double zet = ( rho_up - rho_dn ) / rhotot;
  double g = 0.5 * ( pow(1.0+zet, third2) + pow(1.0-zet, third2) );
  double fk = pi32third * rh13;
  double rs = alpha / fk;
  double twoksg = 2.0 * sqrt( four_over_pi * fk ) *g;
  double t = grad / ( twoksg * rhotot );

  double rtrs = sqrt(rs);
  gPW92 ( 0.0310907, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294,
          rtrs, &eu, &eurs );
  gPW92 ( 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517,
          rtrs, &ep, &eprs );
  gPW92 ( 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671,
          rtrs, &alfm, &alfrsm );
  double z4 = zet * zet * zet * zet;
  double f = (pow(1.0+zet,third4)+pow(1.0-zet,third4)-2.0)/gam;
  ec = eu * ( 1.0 - f * z4 ) + ep * f * z4 - alfm * f * (1.0-z4) / fzz;

  /* LSD potential from [c] (A1) */
  /* ecrs = d ec / d rs [c] (A2) */
  double ecrs = eurs * ( 1.0 - f * z4 ) + eprs * f * z4
                - alfrsm * f * (1.0-z4)/fzz;
  double fz = third4 * ( pow(1.0+zet,third) - pow(1.0-zet,third))/gam;
  double eczet = 4.0 * (zet*zet*zet) * f * ( ep - eu + alfm/fzz ) +
          fz * ( z4 * ep - z4 * eu - (1.0-z4) * alfm/fzz );
  double comm = ec - rs * ecrs * third - zet * eczet;
  vc1_up = comm + eczet;
  vc1_dn = comm - eczet;

  /* SCAN correlation energy */
  /* b = A of [a] (8) */

  double g3 = g * g * g;
  double pon = - ec / (g3 * gamma);
  double b = delt / ( exp ( pon ) - 1.0 );
  double b2 = b * b;
  double t2 = t * t;
  double t4 = t2 * t2;
  double q4 = 1.0 + b * t2;
  double q5 = q4 + b2 * t4;
  double h = g3 * gamma * log ( 1.0 + delt * q4 * t2 / q5 );

  /* Energy done, now the potential, using appendix E of [b] */

  double g4 = g3 * g;
  double t6 = t4 * t2;
  double rsthrd = rs * third;
  double gz = ( pow ( (1.0+zet)*(1.0+zet) + eta, sixthm ) -
         pow ( (1.0-zet)*(1.0-zet) + eta, sixthm ) ) * third;
  double fac = delt / b + 1.0;
  double bg = -3.0 * b2 * ec * fac / ( bet * g4 );
  double bec = b2 * fac / ( bet * g3 );
  double q8 = q5 * q5 + delt * q4 * q5 * t2;
  double q9 = 1.0 + 2.0 * b * t2;
  double hb = - bet * g3 * b * t6 * ( 2.0 + b * t2 ) / q8;
  double hrs = -rsthrd * hb * bec * ecrs;
  double hzed = 3.0 * gz * h / g + hb * ( bg * gz + bec * eczet );
  double ht = 2.0 * bet * g3 * q9 / q8;

  double ccomm = h + hrs - t2 * ht * seven_sixth;
  double pref = hzed - gz * t2 * ht / g;

  ccomm -= pref * zet;

  vc1_up += ccomm + pref;
  vc1_dn += ccomm - pref;
  vc2 = - ht / ( rhotot * twoksg * twoksg );

  *exc_up = x_coeff_ * ex_up + c_coeff_ * ( ec + h );
  *exc_dn = x_coeff_ * ex_dn + c_coeff_ * ( ec + h );
  *vxc1_up = x_coeff_ * vx1_up + c_coeff_ * vc1_up;
  *vxc1_dn = x_coeff_ * vx1_dn + c_coeff_ * vc1_dn;
  *vxc2_upup = x_coeff_ * 2 * vx2_up + c_coeff_ * vc2;
  *vxc2_dndn = x_coeff_ * 2 * vx2_dn + c_coeff_ * vc2;
  *vxc2_updn = c_coeff_ * vc2;
  *vxc2_dnup = c_coeff_ * vc2;
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
