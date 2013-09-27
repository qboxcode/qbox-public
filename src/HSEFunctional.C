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
// Heyd et al.,      J. Chem. Phys. 118, 8207 (2003)
// Heyd, Scuseria    J. Chem. Phys. 120, 7274 (2004)
// Krukau et al.,    J. Chem. Phys. 125, 224106 (2006)
// The PBE exchange hole J is defined here
// Ernzerhof, Perdew J. Chem. Phys. 109, 3313 (1998)
// Parts of this code are taken from the implementation in FLAPW
// Schlipf et al.    Phys. Rev. B 84, 125142 (2011)
//
////////////////////////////////////////////////////////////////////////////////
//
// Author: Martin Schlipf (2013)
// Contact: martin.schlipf@gmail.com
//

#include "HSEFunctional.h"
#include "ExponentialIntegral.h"
#include <cassert>
#include <cmath>
#include <numeric>

using namespace std;

// screening parameter of the HSE functional
const double omega = 0.11;
// parameters of the exchange hole
const double A = 1.0161144, B = -0.37170836, C = -0.077215461, D = 0.57786348,
    E = -0.051955731;

// constructor
HSEFunctional::HSEFunctional(const vector<vector<double> > &rhoe) :
  x_coeff_(0.75), c_coeff_(1.0) {
  // nonmagnetic or magnetic
  _nspin = rhoe.size();
  // in magnetic calculations both density arrays must have the same size
  if (_nspin > 1)
    assert(rhoe[0].size() == rhoe[1].size());

  // allocate arrays
  _np = rhoe[0].size();
  if (_nspin == 1) {
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
  } else {
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

// helper function to calculate the integral
//
// inf
//   /                   /        2                  \
//  |         / w y \   |  A   -Dy           A        |     /        2  2 \
//  | dy Erfc( ----- )  | --- e     - --------------  | Exp( - H(s) s  y   )
//  |         \ k   /   |  y            /    4   2\   |     \             /
// /             F       \            y( 1 + - Ay  ) /
//  0                                   \    9    /
//
// complementary error function is approximated by polynomial
//             8
//           -----             2
//            \        i  - b x
// erfc(x) =   )   a  x  e          for x < 14
//            /     i
//           -----
//           i = 1
// and by exp( -2 x^2 ) above x = 14
// then the integrals with the first part of the exchange hole
// become analytically solvable
void approximateIntegral(const double omega_kF, const double Hs2,
    const double D_term, const double dHs2_ds, double *appInt,
    double *dAppInt_ds, double *dAppInt_dkF) {

  // constant parameterization of error function
  const double a[] = { 1.0, -1.128223946706117, 1.452736265762971,
      -1.243162299390327, 0.971824836115601, -0.568861079687373,
      0.246880514820192, -0.065032363850763, 0.008401793031216 };
  const double b = 1.455915450052607, cutoff = 14.0;

  // helper variables
  const double SQRT_PI = sqrt(M_PI);
  const double A_2 = 0.5 * A, r9_4A = 2.25 / A, sqrtA = sqrt(A);
  const double w2 = omega_kF * omega_kF;
  const double bw2 = b * w2;
  const double r2bw2 = 2.0 * bw2;
  const double bw2_Hs2 = bw2 + Hs2;
  const double bw2_D_term = bw2 + D_term;

  if (bw2_Hs2 < cutoff) {

    // small x limit

    // ### begin calculation of integrals ###

    // inf
    //   /        /        2                  \      /         2              \
    //  |     n  |  A   -Dy           A        |    |    /    w         2\   2 |
    //  | dy y   | --- e     - --------------  | Exp| - (  b --2-- + H s  ) y  |
    //  |        |  y            /    4   2\   |    |    \    k          /     |
    // /          \            y( 1 + - Ay  ) /      \         F              /
    //  0                        \    9    /
    //

    // maximum n in above expression
    // note: to calculate higher derivatives you may need to increase this value
    const int no_integral = 11;

    // Calculate more helper variables
    const double bw2_Hs2_Sqr = bw2_Hs2 * bw2_Hs2;
    const double bw2_Hs2_Cub = bw2_Hs2 * bw2_Hs2_Sqr;
    const double sqrt_bw2_Hs2 = sqrt(bw2_Hs2);
    const double bw2_D_term_Sqr = bw2_D_term * bw2_D_term;
    const double bw2_D_term_Cub = bw2_D_term * bw2_D_term_Sqr;
    const double bw2_D_term_Tet = bw2_D_term_Sqr * bw2_D_term_Sqr;
    const double sqrt_bw2_D_term = sqrt(bw2_D_term);

    // arg = 9/4 * (b w^2/kF^2 + H s^2) / A
    const double arg = r9_4A * bw2_Hs2;
    const double sqrt_arg = sqrt(arg);
    const double r1_arg = 1.0 / arg;

    // calculate e^(arg), E1(arg), and erfc(sqrt(arg))
    const double exp_arg = exp(arg);
    const double exp_erfc = exp_arg * erfc(sqrt_arg);

    // evaluate exponenential integral
    double term2 = (arg < util::series_cutoff) ? exp_arg * util::E1(arg)
        : util::gauss_laguerre(arg);

    // allocate array
    vector<double> integral(no_integral);

    // The n = 0 integral is
    // A/2 ( ln((b (w/kF)^2 + H s^2) / (b (w/kF)^2 + D + H s^2))
    //     + e^(arg) E1(arg) )
    integral[0] = A_2 * (log(bw2_Hs2 / bw2_D_term) + term2);

    // Calculate now all even n's by successive derivation
    // The log(...) term gives term proportional to 1/(b (w/kF)^2 + D + H s^2)^i
    // The e^(arg) E1(arg) term reproduces itself with a prefactor and
    // generates an additional 1/arg term which produces higher 1/arg^i terms
    // when further deriviated
    double term1 = A_2 / bw2_D_term;
    double factor2 = -1.125;
    double arg_n = -1.0 / arg;
    integral[2] = term1 + factor2 * term2;

    for (int i = 1; i < no_integral / 2; i++) {
      term1 = term1 / bw2_D_term * static_cast<double> (i);
      factor2 = -factor2 * r9_4A;
      term2 = term2 + arg_n;

      integral[2 * (i + 1)] = term1 + factor2 * term2;

      arg_n = -arg_n * static_cast<double> (i) / arg;
    }

    // The n = 1 integral is
    // A/2 ( sqrt(pi) / sqrt( b (w/kF)^2 + D + H s2 )
    // - 3/4 sqrt(A) pi e^(arg) erfc(sqrt(arg))
    term1 = A_2 * SQRT_PI / sqrt_bw2_D_term;
    term2 = M_PI * exp_erfc;
    factor2 = -0.75 * sqrtA;

    integral[1] = term1 + factor2 * term2;

    // Calculate now all uneven n's by successive derivation
    // The 1 / sqrt(...) term gives higher orders of 1 / (...)^((2i+1)/2)
    // The e^(arg) erfc(sqrt(arg)) term reproduces itself with a prefactor
    // and generates an additional 1/sqrt(arg) term which produces higher
    // 1/(arg)^((2i+1)/2) terms when further deriviated
    double sum_term = -1.125 * SQRT_PI / sqrt_bw2_Hs2;
    double add_term = sum_term;
    double half_i2_1 = -0.5;

    for (int i = 3; i < no_integral; i += 2) {
      factor2 = -factor2 * r9_4A;
      term1 = -term1 * half_i2_1 / bw2_D_term;
      integral[i] = term1 + term2 * factor2 + sum_term;

      add_term = -add_term * half_i2_1 / bw2_Hs2;
      sum_term = -sum_term * r9_4A + add_term;
      half_i2_1 = half_i2_1 - 1.0;
    }

    // ### end calculation of integrals ###

    const int no_coeff = sizeof(a) / sizeof(double);

    // allocate vector
    vector<double> a_wi, ai_wi;
    // will contain a[i] * (w/kF)^i
    a_wi.reserve(no_coeff);
    // will contain a[i] * i * (w/kF)^i
    ai_wi.reserve(no_coeff);

    // initialize array
    double wi = 1.0;
    for (int i = 0; i < no_coeff; i++) {
      a_wi.push_back(a[i] * wi);
      ai_wi.push_back(a_wi[i] * i);
      wi *= omega_kF;
    }

    //  combine the solutions of the integrals with the appropriate prefactors
    *appInt = inner_product(a_wi.begin(), a_wi.end(), integral.begin(), 0.0);
    // for derivative shift integral index by 2
    const double dotpr = inner_product(a_wi.begin(), a_wi.end(),
        integral.begin() + 2, 0.0);
    *dAppInt_ds = -dotpr * dHs2_ds;
    *dAppInt_dkF = -inner_product(ai_wi.begin(), ai_wi.end(), integral.begin(),
        r2bw2 * dotpr);

  } else {

    // large x limit

    // inf
    //   /    /        2                  \      /         2              \
    //  |    |  A   -Dy           A        |    |    /    w         2\   2 |
    //  | dy | --- e     - --------------  | Exp| - (  2 --2-- + H s  ) y  |
    //  |    |  y            /    4   2\   |    |    \    k          /     |
    // /      \            y( 1 + - Ay  ) /      \         F              /
    //  0                    \    9    /

    const double r2w2 = 2.0 * w2;
    const double r4w2 = 4.0 * w2;
    const double r2w2_Hs2 = r2w2 + Hs2;
    const double r2w2_D_term = r2w2 + D_term;
    const double arg = r9_4A * r2w2_Hs2;
    const double exp_e1 = util::gauss_laguerre(arg);
    *appInt = A_2 * (log(r2w2_Hs2 / r2w2_D_term) + exp_e1);
    const double dAppInt_dh = -A_2 / r2w2_D_term + 1.125 * exp_e1;
    *dAppInt_ds = dAppInt_dh * dHs2_ds;
    *dAppInt_dkF = -dAppInt_dh * r4w2;

  }

}

// helper function to calculate the HSE enhancement factor
// and its derivative w.r.t. s and k_F
//                inf
//                  /
//  HSE         8  |                  w y
// F (s,k ) = - -  | dy J(s,y) erfc( ----- )
//  x    F      9  |                   k
//                /                     F
//                 0
// where s:   reduced density gradient
//       k_F: Fermi vector
//       w:   screening of the functional
// the PBE exchange hole is defined by the following expression
//           /
//          |    A          1          /  A                  2         2
// J(s,y) = | - --2- ------------2- + |  --2- + B + C [ 1 + s  F(s) ] y
//          |    y    1 + (4/9)Ay      \  y
//           \
//                                           2 \      2       2
//                       2         4 \   -D y   |  - s  H(s) y
//            + E [ 1 + s  G(s) ] y   | e       | e
//                                   /          |
//                                             /
// see references for a definition of the functions F(s), G(s), and H(s)
void HSE_enhance(const double s, const double kF, const double w, double *fx,
    double *dfx_ds, double* dfx_dkf) {

  // prefactor for exchange hole
  const double r8_9 = 8.0 / 9.0;

  // derived quantities
  const double SQRT_PI = sqrt(M_PI);
  const double sqrtA = sqrt(A), r9_4A = 2.25 * A;

  // evaluate
  //                 2       4
  //             a1 s  + a2 s
  // H(s) = ---------4------5------6-
  //         1 + a3 s + a4 s + a5 s
  // and
  //     2                 3         5            3        4        5
  // d( s H(s) )   ( 4 a1 s  + 6 a2 s ) - ( 4 a3 s + 5 a4 s + 6 a5 s ) * H(s)
  // ----------- = ---------------------4------5------6----------------------
  //     d s                    1 + a3 s + a4 s + a5 s

  const double a1 = 0.00979681, a2 = 0.0410834, a3 = 0.187440, a4 = 0.00120824,
      a5 = 0.0347188;

  // helper variables
  const double s2 = s * s;
  const double s3 = s2 * s;
  const double s4 = s2 * s2;
  const double s5 = s3 * s2;
  const double s6 = s3 * s3;

  // calculate numerator and reciprocal of denominator
  const double numerator = a1 * s2 + a2 * s4;
  const double r_denom = 1.0 / (1.0 + a3 * s4 + a4 * s5 + a5 * s6);
  // helper for derivatives
  const double first = 4.0 * a1 * s3 + 6.0 * a2 * s5;
  const double second = 4.0 * a3 * s3 + 5.0 * a4 * s4 + 6.0 * a5 * s5;

  // put everything together
  const double H = numerator * r_denom;
  const double Hs2 = H * s2;
  const double dHs2_ds = (first - second * H) * r_denom;

  // evaluate
  // F(s) = 6.475 H(s) + 0.4797
  // d( s^2 F(s) )         d( s^2 H(s) )
  // ------------- = 6.475 ------------- + 2 * 0.4797
  //     d s
  const double slope = 6.475;
  const double intercept = 0.4797;

  const double F = slope * H + intercept;
  const double Fs2 = F * s2;
  const double dFs2_ds = slope * dHs2_ds + 2 * intercept;

  // evaluate alpha and beta
  // alpha = part1 - part2
  //               __  2
  //           15 VPi s
  // beta = -----------------
  //                  2  7/2
  //         16 (D + H s  )
  // with
  //
  //                          2      2           2 2          2 3
  //          __ 15E + 6C(1+Fs )(D+Hs ) + 4B(D+Hs )  + 8A(D+Hs )
  // part1 = VPi ------------------------------------------------
  //                                    2  7/2
  //                        16 ( D + H s  )
  // and
  //                                            ________
  //                       /     2 \       /   /     2   \
  //         3 Pi  ___    | 9 H s   |     |   / 9 H s     |
  // part2 = ---- V A  Exp| ------- | Erfc|  /  -------   |
  //          4           |   4 A   |     | V     4 A     |
  //                       \       /       \             /
  // and the derivatives with respect to s

  // calculate the helper variables
  const double AHs2_1_2 = sqrtA * sqrt(Hs2);
  const double AHs2_3_2 = AHs2_1_2 * A * Hs2;
  const double r1_Fs2 = 1.0 + Fs2;
  const double D_Hs2 = D + Hs2;
  const double D_Hs2Sqr = D_Hs2 * D_Hs2;
  const double D_Hs2Cub = D_Hs2 * D_Hs2Sqr;
  const double D_Hs2_5_2 = D_Hs2Sqr * sqrt(D_Hs2);
  const double D_Hs2_7_2 = D_Hs2_5_2 * D_Hs2;
  const double D_Hs2_9_2 = D_Hs2_5_2 * D_Hs2Sqr;
  const double D_Hs2_11_2 = D_Hs2_5_2 * D_Hs2Cub;

  // part 1 and derivatives w.r.t. Hs^2 and Fs^2
  const double part1 = SQRT_PI * (15.0 * E + 6.0 * C * r1_Fs2 * D_Hs2 + 4.0 * B
      * D_Hs2Sqr + 8.0 * A * D_Hs2Cub) / (16.0 * D_Hs2_7_2);
  const double dpart1_dh = -SQRT_PI * (105.0 * E + 30.0 * C * r1_Fs2 * D_Hs2
      + 12.0 * B * D_Hs2Sqr + 8.0 * A * D_Hs2Cub) / (32.0 * D_Hs2_9_2);
  const double dpart1_df = SQRT_PI * 0.375 * C / D_Hs2_5_2;

  // part 2 and derivative w.r.t. Hs^2
  const double arg1 = r9_4A * Hs2;
  const double arg2 = sqrt(arg1);
  const double exp_erfc = exp(arg1) * erfc(arg2);
  const double part2 = 0.75 * M_PI * sqrtA * exp_erfc;
  const double dpart2_dh = 0.75 * M_PI * sqrtA * (r9_4A * exp_erfc - 1.5
      / (SQRT_PI * AHs2_1_2));

  // combine parts and derivatives
  const double alpha = part1 - part2;
  const double dalpha_dh = dpart1_dh - dpart2_dh;
  const double dalpha_df = dpart1_df;

  // calculate beta / s^2, its derivative w.r.t. Hs^2 and E * beta
  const double beta_s2 = 0.9375 * SQRT_PI / D_Hs2_7_2;
  const double Ebeta = E * beta_s2 * s2;
  const double dbeta_dh = -3.5 * beta_s2 / D_Hs2;

  // combine alpha and beta to function G
  //       3 Pi / 4 + alpha
  // G = - ----------------
  //             beta
  const double r3Pi_4_alpha = 0.75 * M_PI + alpha;
  const double G = -r3Pi_4_alpha / Ebeta;

  // calculate derivative w.r.t. s
  //
  //             /  /3 Pi    \  d(b/s^2)   /           d a    \  d(Hs^2)       d a   d(Fs^2)
  //            <  ( ---- + a ) -------   /  b/s^2  - -------  > -------  -  ------- -------
  // d (Gs^2)    \  \ 4      /  d(Hs^2)  /            d(Hs^2) /    d s       d(Fs^2)   ds
  // -------- = ----------------------------------------------------------------------------
  //    ds                                     E b/s^2
  //
  // notice that alpha and beta are abbreviated by a and b in this equation
  const double Ebeta_s2 = E * beta_s2;
  const double dGs2_ds = ((r3Pi_4_alpha * dbeta_dh / beta_s2 - dalpha_dh)
      * dHs2_ds - dalpha_df * dFs2_ds) / Ebeta_s2;

  // helper variables for the integration of the exchange hole
  const double C_term = C * (1 + s2 * F);
  const double dCt_ds = C * dFs2_ds;
  const double E_term = E * (1 + s2 * G);
  const double dEt_ds = E * dGs2_ds;
  const double D_term = D + Hs2;
  const double r1_D_term = 1.0 / D_term;
  const double r1_kF = 1.0 / kF;
  const double w_kF = w * r1_kF;
  const double w_kF_Sqr = w_kF * w_kF;
  const double r1_a = w_kF_Sqr / D_term;
  const double r1_sqrta = sqrt(r1_a);

  // approximate the integral using an expansion of the error function
  double appInt, dAppInt_ds, dAppInt_dkF;
  approximateIntegral(w_kF, Hs2, D_term, dHs2_ds, &appInt, &dAppInt_ds,
      &dAppInt_dkF);

  // Calculate the integrals
  //
  // inf
  //   /                       2   2      /     \
  //  |     2 n + 1   -(D+H(s)s ) y      |  w    |
  //  | dy y         e               Erfc| --- y |
  //  |                                  |  k    |
  // /                                    \  F  /
  //  0
  // we use that
  // inf                                                      n
  //   /                  2                           /     -----                                  \
  //  |     2 n + 1  - a y                 n   n!    |       \     (2n - 1)!!             -(2n + 1) |
  //  | dy y        e      Erfc(y) = ( -1 )  ----n-  |  1 -   )    ---------- sqrt( a + 1 )         |
  //  |                                       2 a    |       /       (2n)!!                         |
  // /                                                \     -----                                  /
  //  0                                                     m = 0
  // with a = (D + H(s)s^2) k_F^2 / w^2
  // and n!! = n * (n-2) * (n-4) ...; if (n <= 0) := 1

  // allocate array
  vector<double> intYExpErfc(0);
  intYExpErfc.reserve(4);

  // initialize
  double prefactor = 0.5 / D_term; // note w_kF cancels in transformation
  double summand = r1_sqrta;
  double sum = 1 - summand;
  intYExpErfc.push_back(prefactor * sum);

  // calculate higher n integrals
  for (int i = 1; i < 4; i++) {

    // update values
    prefactor *= -static_cast<double> (i) * r1_D_term;
    summand *= static_cast<double> (2 * i - 1) / static_cast<double> (2 * i)
        * r1_a;
    sum -= summand;

    intYExpErfc.push_back(prefactor * sum);

  }

  // Calculate the integrals
  //
  //       inf
  //         /                 /                    2      \
  //   2    |     2 (n+1)     |          2   2     w     2  |
  // -----  | dy y         Exp| -(D+H(s)s ) y  - ------ y   |
  //  ____  |                 |                   k ^2      |
  // V Pi  /                   \                   F       /
  //        0
  // for n = 0, 1, 2
  //
  const double r1_arg = 1.0 / (D_term + w_kF_Sqr);
  // allocate array
  vector<double> intYGauss(3);
  intYGauss[0] = 0.5 * sqrt(r1_arg) * r1_arg;
  intYGauss[1] = 1.5 * intYGauss[0] * r1_arg;
  intYGauss[2] = 2.5 * intYGauss[1] * r1_arg;

  // put everything together

  // Calculate the integral
  //  inf
  //    /                 /      \
  //   |                 |  w     |
  //   | dy y J(s,y) Erfc| ---- y |
  //   |                 |  k     |
  //  /                   \  F   /
  //   0
  // where J(s, y) is the exchange hole defined in the references
  // the exchange factor is proportional to this integral
  *fx = -r8_9 * (appInt + B * intYExpErfc[0] + C_term * intYExpErfc[1] + E_term
      * intYExpErfc[2]);

  // Calculate the derivatives with respect to s using that the derivatative of the integral
  // yields higher orders of the same kind of integral intY1 -> -intY3 -> intY5 ... times
  // the derivative of the exponent
  *dfx_ds = -r8_9 * (dAppInt_ds - (B * intYExpErfc[1] + C_term * intYExpErfc[2]
      + E_term * intYExpErfc[3]) * dHs2_ds + dCt_ds * intYExpErfc[1] + dEt_ds
      * intYExpErfc[2]);
  *dfx_dkf = -r8_9 * r1_kF * (w_kF * (B * intYGauss[0] + C_term * intYGauss[1]
      + E_term * intYGauss[2]) + dAppInt_dkF);

}

void HSEFunctional::setxc(void) {
  // dummy
}
