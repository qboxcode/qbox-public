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
// RSHFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Range-separated hybrid functional (RSH)
// J.H.Skone et al. Phys. Rev. B93, 235106 (2016)
// RSH is defined by alpha_RSH, beta_RSH, mu_RSH
// sigma = alpha_RSH * rho(r,r') * erf(r-r')/(r-r') +
//         beta_RSH * rho(r,r') * erfc(r-r')/(r-r') +
//         (1 - alpha_RSH) * Vx_LR(r,mu_RSH) +
//         (1 - beta_RSH) * Vx_SR(r,mu_RSH)
// The HSE functional is obtained using alpha_RSH=0, beta_RSH=0.25, mu_RSH=0.11
// Heyd et al.,      J. Chem. Phys. 118, 8207 (2003)
// Heyd, Scuseria    J. Chem. Phys. 120, 7274 (2004)
// Krukau et al.,    J. Chem. Phys. 125, 224106 (2006)
// The PBE exchange hole J is defined here
// Ernzerhof, Perdew J. Chem. Phys. 109, 3313 (1998)
// Parts of this code are taken from the implementation in FLAPW
// Schlipf et al.    Phys. Rev. B 84, 125142 (2011)
//
////////////////////////////////////////////////////////////////////////////////

#include "RSHFunctional.h"
#include "ExponentialIntegral.h"
#include <cassert>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace std;

// parameters of the exchange hole
const double A = 1.0161144, B = -0.37170836, C = -0.077215461, D = 0.57786348,
  E = -0.051955731;

// constructor
RSHFunctional::RSHFunctional(const vector<vector<double> > &rhoe,
  double alpha_RSH, double beta_RSH, double mu_RSH):
  alpha_RSH_(alpha_RSH), beta_RSH_(beta_RSH), mu_RSH_(mu_RSH), omega(mu_RSH),
  x_coeff_(1.0-beta_RSH), c_coeff_(1.0)
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
void RSHFunctional::approximateIntegral(const double omega_kF, const double Hs2,
  const double D_term, const double dHs2_ds, double *appInt,
  double *dAppInt_ds, double *dAppInt_dkF)
{

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

  if ( bw2_Hs2 < cutoff )
  {

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
    const double sqrt_bw2_Hs2 = sqrt(bw2_Hs2);
    const double sqrt_bw2_D_term = sqrt(bw2_D_term);

    // arg = 9/4 * (b w^2/kF^2 + H s^2) / A
    const double arg = r9_4A * bw2_Hs2;
    const double sqrt_arg = sqrt(arg);

    // calculate e^(arg), E1(arg), and erfc(sqrt(arg))
    const double exp_arg = exp(arg);
    const double exp_erfc = exp_arg * erfc(sqrt_arg);

    // evaluate exponenential integral
    double term2 = ( arg < util::series_cutoff ) ? exp_arg * util::E1(arg)
      : util::gauss_laguerre(arg);

    // allocate array
    vector<double> integral(no_integral);

    // The n = 0 integral is
    // A/2 ( ln((b (w/kF)^2 + H s^2) / (b (w/kF)^2 + D + H s^2))
    //     + e^(arg) E1(arg) )
    integral[0] = A_2 * ( log(bw2_Hs2 / bw2_D_term) + term2 );

    // Calculate now all even n's by successive derivation
    // The log(...) term gives term proportional to 1/(b (w/kF)^2 + D + H s^2)^i
    // The e^(arg) E1(arg) term reproduces itself with a prefactor and
    // generates an additional 1/arg term which produces higher 1/arg^i terms
    // when further deriviated
    double term1 = A_2 / bw2_D_term;
    double factor2 = -1.125;
    double arg_n = -1.0 / arg;
    integral[2] = term1 + factor2 * term2;

    for ( int i = 1; i < no_integral / 2; i++ )
    {
      term1 = i * term1 / bw2_D_term;
      factor2 = -factor2 * r9_4A;
      term2 = term2 + arg_n;

      integral[2 * ( i + 1 )] = term1 + factor2 * term2;

      arg_n = -arg_n * i / arg;
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

    for ( int i = 3; i < no_integral; i += 2 )
    {
      factor2 = -factor2 * r9_4A;
      term1 = -term1 * half_i2_1 / bw2_D_term;
      integral[i] = term1 + term2 * factor2 + sum_term;

      add_term = -add_term * half_i2_1 / bw2_Hs2;
      sum_term = -sum_term * r9_4A + add_term;
      half_i2_1 = half_i2_1 - 1.0;
    }

    // ### end calculation of integrals ###

    const int no_coeff = sizeof( a ) / sizeof(double);

    // allocate vector
    vector<double> a_wi, ai_wi;
    // will contain a[i] * (w/kF)^i
    a_wi.reserve(no_coeff);
    // will contain a[i] * i * (w/kF)^i
    ai_wi.reserve(no_coeff);

    // initialize array
    double wi = 1.0;
    for ( int i = 0; i < no_coeff; i++ )
    {
      a_wi.push_back(a[i] * wi);
      ai_wi.push_back(a_wi[i] * i);
      wi *= omega_kF;
    }

    //  combine the solutions of the integrals with the appropriate prefactors
    *appInt = inner_product(a_wi.begin(),a_wi.end(),integral.begin(),0.0);
    // for derivative shift integral index by 2
    const double dotpr = inner_product(a_wi.begin(),a_wi.end(),integral.begin()
      + 2,0.0);
    *dAppInt_ds = -dotpr * dHs2_ds;
    *dAppInt_dkF = -inner_product(ai_wi.begin(),ai_wi.end(),integral.begin(),
      -r2bw2 * dotpr);

  }
  else
  {

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
    *appInt = A_2 * ( log(r2w2_Hs2 / r2w2_D_term) + exp_e1 );
    const double dAppInt_dh = -A_2 / r2w2_D_term + 1.125 * exp_e1;
    *dAppInt_ds = dAppInt_dh * dHs2_ds;
    *dAppInt_dkF = -dAppInt_dh * r4w2;

  }

}

// helper function to calculate the HSE enhancement factor
// and its derivative w.r.t. s and k_F
//                inf
//                  /
//  HSE         8  |                    w y
// F (s,k ) = - -  | dy y J(s,y) erfc( ----- )
//  x    F      9  |                    k
//                /                      F
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
void RSHFunctional::RSH_enhance(const double s_inp, const double kF,
  const double w, double *fx, double *dfx_ds, double* dfx_dkf)
{
  // Correction of the reduced gradient to ensure Lieb-Oxford bound
  // If a large value of s would violate the Lieb-Oxford bound, the
  // value of s is reduced, so that this condition is fullfilled
  const double s_thresh = 8.3, s_max = 8.572844, s_chg = 18.796223;
  const bool correction = s_inp > s_thresh;
  const double s = ( correction ) ? s_max - s_chg / ( s_inp * s_inp ) : s_inp;

  // sanity check
  assert(s > 0);

  // prefactor for exchange hole
  const double r8_9 = 8.0 / 9.0;

  // derived quantities
  const double SQRT_PI = sqrt(M_PI);
  const double sqrtA = sqrt(A);
  const double r9_4A = 2.25 / A;

  // evaluate
  //                 2       4
  //             a1 s  + a2 s
  // H(s) = ---------4------5------6-
  //         1 + a3 s + a4 s + a5 s
  // and
  //   2                 3         5            3        4        5          2
  // d(s H(s))   ( 4 a1 s  + 6 a2 s ) - ( 4 a3 s + 5 a4 s + 6 a5 s ) * H(s) s
  // --------- = ---------------------4------5------6--------------------------
  //    d s                    1 + a3 s + a4 s + a5 s

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
  const double r_denom = 1.0 / ( 1.0 + a3 * s4 + a4 * s5 + a5 * s6 );
  // helper for derivatives
  const double first = 4.0 * a1 * s3 + 6.0 * a2 * s5;
  const double second = 4.0 * a3 * s3 + 5.0 * a4 * s4 + 6.0 * a5 * s5;

  // put everything together
  const double H = numerator * r_denom;
  const double Hs2 = H * s2;
  const double dHs2_ds = ( first - second * Hs2 ) * r_denom;

  // evaluate
  // F(s) = Int + Sl * H(s)
  // d( s^2 F(s) )                     d( s^2 H(s) )
  // ------------- = 2 * Int *s + Sl * -------------
  //     d s                                d s
  // Int - intercept
  // Sl  - slope
  const double slope = 6.4753871;
  const double intercept = 0.4796583;

  const double F = slope * H + intercept;
  const double Fs2 = F * s2;
  const double dFs2_ds = slope * dHs2_ds + 2 * intercept * s;

  // evaluate alpha and beta
  // alpha = part1 - part2
  //               __  2
  //           15 VPi s
  // beta = -----------------
  //                    2  7/2
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
  const double r1_Fs2 = 1.0 + Fs2;
  const double D_Hs2 = D + Hs2;
  const double D_Hs2Sqr = D_Hs2 * D_Hs2;
  const double D_Hs2Cub = D_Hs2 * D_Hs2Sqr;
  const double D_Hs2_5_2 = D_Hs2Sqr * sqrt(D_Hs2);
  const double D_Hs2_7_2 = D_Hs2_5_2 * D_Hs2;
  const double D_Hs2_9_2 = D_Hs2_5_2 * D_Hs2Sqr;

  // part 1 and derivatives w.r.t. Hs^2 and Fs^2
  const double part1 = SQRT_PI * ( 15.0 * E + 6.0 * C * r1_Fs2 * D_Hs2 + 4.0
    * B * D_Hs2Sqr + 8.0 * A * D_Hs2Cub ) / ( 16.0 * D_Hs2_7_2 );
  const double dpart1_dh = -SQRT_PI * ( 105.0 * E + 30.0 * C * r1_Fs2 * D_Hs2
    + 12.0 * B * D_Hs2Sqr + 8.0 * A * D_Hs2Cub ) / ( 32.0 * D_Hs2_9_2 );
  const double dpart1_df = SQRT_PI * 0.375 * C / D_Hs2_5_2;

  // part 2 and derivative w.r.t. Hs^2
  const double arg1 = r9_4A * Hs2;
  const double arg2 = sqrt(arg1);
  const double exp_erfc = exp(arg1) * erfc(arg2);
  const double part2 = 0.75 * M_PI * sqrtA * exp_erfc;
  const double dpart2_dh = 0.75 * M_PI * sqrtA * ( r9_4A * exp_erfc - 1.5
    / ( SQRT_PI * AHs2_1_2 ) );

  // combine parts and derivatives
  const double alpha = part1 - part2;
  const double dalpha_dh = dpart1_dh - dpart2_dh;
  const double dalpha_df = dpart1_df;

  // calculate beta / s^2, its derivative w.r.t. Hs^2 and E * beta
  const double Ebeta_s2 = E * 0.9375 * SQRT_PI / D_Hs2_7_2;
  const double dEbeta_dh = -3.5 * Ebeta_s2 / D_Hs2;

  // combine alpha and beta to function G
  //       3 Pi / 4 + alpha
  // G = - ----------------
  //            E beta
  const double r3Pi_4_alpha = 0.75 * M_PI + alpha;
  // Gs2 = G * s^2
  const double Gs2 = -r3Pi_4_alpha / Ebeta_s2;

  // calculate derivative w.r.t. s
  //
  //             /  /3 Pi    \  d(b/s^2)   /           d a    \  d(Hs^2)
  //            <  ( ---- + a ) -------   /  b/s^2  - -------  > -------
  // d (Gs^2)    \  \ 4      /  d(Hs^2)  /            d(Hs^2) /    d s
  // -------- = --------------------------------------------------------
  //    ds                             E b/s^2
  //
  //
  //        d a   d(Fs^2)
  //      ------- -------
  //      d(Fs^2)   ds
  //  - -----------------
  //        E b/s^2
  //
  // notice that alpha and beta are abbreviated by a and b in this equation
  const double dGs2_ds = ( ( r3Pi_4_alpha * dEbeta_dh / Ebeta_s2 - dalpha_dh ) *
    dHs2_ds - dalpha_df * dFs2_ds ) / Ebeta_s2;

  // helper variables for the integration of the exchange hole
  const double C_term = C * ( 1 + s2 * F );
  const double dCt_ds = C * dFs2_ds;
  const double E_term = E * ( 1 + Gs2 );
  const double dEt_ds = E * dGs2_ds;
  const double D_term = D + Hs2;
  const double r1_D_term = 1.0 / D_term;
  const double r1_kF = 1.0 / kF;
  const double w_kF = w * r1_kF;
  const double w_kF_Sqr = w_kF * w_kF;

  // approximate the integral using an expansion of the error function
  double appInt, dAppInt_ds, dAppInt_dkF;
  approximateIntegral(w_kF,Hs2,D_term,dHs2_ds,&appInt,&dAppInt_ds,&dAppInt_dkF);

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
  // inf
  //   /                  2
  //  |     2 n + 1  - q y              n!
  //  | dy y        e      Erfc(y) = ----n+1-  *
  //  |                               2 q
  // /
  //  0
  //            n
  //    /     -----                                     \
  //   |       \     (2m - 1)!!  m             -(2m + 1) |
  //   |  1 -   )    ---------- q  sqrt( q + 1 )         |
  //   |       /       (2m)!!                            |
  //    \     -----                                     /
  //          m = 0
  // with q = (D + H(s)s^2) k_F^2 / w^2
  // and n!! = n * (n-2) * (n-4) ...; if (n <= 0) := 1

  // allocate array
  vector<double> intYExpErfc(0);
  intYExpErfc.reserve(4);

  // helper variables
  const double q = D_term / w_kF_Sqr;
  const double q_q_1 = q / ( q + 1 );
  const double sqrtq_1 = sqrt(q + 1);

  // initialize
  double prefactor = 0.5 / D_term; // note w_kF cancels in transformation
  double summand = 1.0 / sqrtq_1;
  double sum = 1 - summand;
  intYExpErfc.push_back(prefactor * sum);

  // calculate higher n integrals
  for ( int i = 1; i < 4; i++ )
  {

    // update values
    prefactor *= i * r1_D_term;
    summand *= (2.0 * i - 1.0) / (2.0 * i) * q_q_1;
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
  const double r1_arg = 1.0 / ( D_term + w_kF_Sqr );
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
  *fx = -r8_9 * ( appInt + B * intYExpErfc[0] + C_term * intYExpErfc[1]
    + E_term * intYExpErfc[2] );

  // Calculate the derivatives with respect to s using that the derivatative
  // of the integral yields higher orders of the same kind of integral
  //  intY1 -> -intY3 -> intY5 ... times the derivative of the exponent
  *dfx_ds = -r8_9 * ( dAppInt_ds - ( B * intYExpErfc[1] + C_term
    * intYExpErfc[2] + E_term * intYExpErfc[3] ) * dHs2_ds + dCt_ds
    * intYExpErfc[1] + dEt_ds * intYExpErfc[2] );
  *dfx_dkf = -r8_9 * r1_kF * ( w_kF * ( B * intYGauss[0] + C_term
    * intYGauss[1] + E_term * intYGauss[2] ) + dAppInt_dkF );

  // if the value of s has been corrected to satisfy Lieb-Oxford bound,
  // derivative must be adjusted as well
  if ( correction )
  {
    *dfx_ds *= 2.0 * s_chg * pow(s_inp,-3);
  }

}

// exchange helper function
// input:
// rho  - charge density
// grad - absolute value of gradient
// a_ex - amount of HF exchange
// w    - screening
// output:
// ex       - exchange energy
// vx1, vx2 - exchange potential such that vx = vx1 + div( vx2 * grad(n) )
void RSHFunctional::RSH_exchange(const double rho, const double grad,
  const double a_ex, const double w, double *ex, double *vx1, double *vx2)
{

  // constants employed in the PBE/HSE exchange
  const double third = 1.0 / 3.0;
  const double third4 = 4.0 / 3.0;
  const double ax = -0.7385587663820224058; /* -0.75*pow(3.0/pi,third) */
  const double um = 0.2195149727645171;
  const double uk = 0.804;
  const double ul = um / uk;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */

  // intialize
  *ex = 0;
  *vx1 = 0;
  *vx2 = 0;

  // very small densities do not contribute
  if ( rho < 1e-18 ) return;

  // LDA exchange energy
  const double rho13 = pow(rho,third);
  const double exLDA = ax * rho13;

  // Fermi wave vector kF = ( 3 * pi^2 n )^(1/3)
  const double kF = pi32third * rho13;
  // reduced density gradient
  const double s = grad / ( 2.0 * kF * rho );

  // calculate PBE enhancement factor
  const double s2 = s * s;
  const double p0 = 1.0 + ul * s2;
  const double fxpbe = 1.0 + uk - uk / p0;
  // fs = (1/s) * d Fx / d s
  const double fs = 2.0 * uk * ul / ( p0 * p0 );
  // calculate HSE enhancement factor and derivatives w.r.t. s and kF
  double fxhse, dfx_ds, dfx_dkf;
  RSH_enhance(s,kF,w,&fxhse,&dfx_ds,&dfx_dkf);

  // calculate exchange energy
  *ex = exLDA * ( ( 1.0 - alpha_RSH_ ) * fxpbe +
                  ( alpha_RSH_ - beta_RSH_ ) * fxhse );

  // calculate potential
  *vx1 = third4 * exLDA * ( fxpbe - s2 * fs - a_ex * ( fxhse - s * dfx_ds
    + 0.25 * kF * dfx_dkf ) );
  *vx2 = -exLDA * ( fs / ( rho * 4.0 * kF * kF ) - a_ex * dfx_ds / ( 2.0 * kF
    * grad ) );

}

////////////////////////////////////////////////////////////////////////////////
//
//  gcor2.c: Interpolate LSD correlation energy
//  as given by (10) of Perdew & Wang, Phys Rev B45 13244 (1992)
//  Translated into C by F.Gygi, Dec 9, 1996
//
////////////////////////////////////////////////////////////////////////////////

void RSHFunctional::gcor2(double a, double a1, double b1, double b2,
  double b3, double b4, double rtrs, double *gg, double *ggrs)
{
  double q0, q1, q2, q3;
  q0 = -2.0 * a * ( 1.0 + a1 * rtrs * rtrs );
  q1 = 2.0 * a * rtrs * ( b1 + rtrs * ( b2 + rtrs * ( b3 + rtrs * b4 ) ) );
  q2 = log(1.0 + 1.0 / q1);
  *gg = q0 * q2;
  q3 = a * ( b1 / rtrs + 2.0 * b2 + rtrs * ( 3.0 * b3 + 4.0 * b4 * rtrs ) );
  *ggrs = -2.0 * a * a1 * q2 - q0 * q3 / ( q1 * ( 1.0 + q1 ) );
}

////////////////////////////////////////////////////////////////////////////////
//
//  calculate correlation energy of the PBE functional
//  K.Burke's modification of PW91 codes, May 14, 1996.
//  Modified again by K.Burke, June 29, 1996, with simpler Fx(s)
//  Translated into C and modified by F.Gygi, Dec 9, 1996.
//
//  input:
//    rho:  density
//    grad: abs(grad(rho))
//  output:
//    exc: exchange-correlation energy per electron
//    vxc1, vxc2 : quantities such that the total exchange potential is:
//
//      vxc = vxc1 + div ( vxc2 * grad(n) )
//
//  References:
//  [a] J.P.Perdew, K.Burke, and M.Ernzerhof,
//      "Generalized gradient approximation made simple,
//      Phys.Rev.Lett. 77, 3865, (1996).
//  [b] J.P.Perdew and Y.Wang, Phys.Rev. B33, 8800 (1986),
//      Phys.Rev. B40, 3399 (1989) (E).
//
////////////////////////////////////////////////////////////////////////////////

void RSHFunctional::PBE_correlation(const double rho, const double grad,
  double *ec, double *vc1, double *vc2)
{
  *ec = 0.0;
  *vc1 = 0.0;
  *vc2 = 0.0;

  if ( rho < 1.e-18  )
  {
    return;
  }

  const double third = 1.0 / 3.0;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
  const double alpha = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
  const double seven_sixth = 7.0 / 6.0;
  const double four_over_pi = 1.27323954473516;
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double bet = 0.06672455060314922; /* see [a] (4) */
  const double delt = bet / gamma;

  // Fermi wave vector kF = ( 3 * pi^2 n )^(1/3)
  const double rho13 = pow(rho,third);
  const double fk = pi32third * rho13;

  /* Find LSD contributions, using [c] (10) and Table I of [c]. */
  /* ec = unpolarized LSD correlation energy */
  /* ecrs = d ec / d rs */
  /* construct ec, using [c] (8) */

  const double rs = alpha / fk;
  const double twoks = 2.0 * sqrt(four_over_pi * fk);
  const double t = grad / ( twoks * rho );

  const double rtrs = sqrt(rs);
  double ecrs;
  gcor2(0.0310907,0.2137,7.5957,3.5876,1.6382,0.49294,rtrs,ec,&ecrs);

  /* LSD potential from [c] (A1) */
  /* ecrs = d ec / d rs [c] (A2) */

  const double vc = *ec - rs * ecrs * third;

  /* PBE correlation energy */
  /* b = A of [a] (8) */

  const double pon = -*ec / gamma;
  const double b = delt / ( exp(pon) - 1.0 );
  const double b2 = b * b;
  const double t2 = t * t;
  const double t4 = t2 * t2;
  const double q4 = 1.0 + b * t2;
  const double q5 = q4 + b2 * t4;
  const double h = gamma * log(1.0 + delt * q4 * t2 / q5);

  // Energy done, now the potential, using appendix E of [b]

  const double t6 = t4 * t2;
  const double rsthrd = rs * third;
  const double fac = delt / b + 1.0;
  const double bec = b2 * fac / bet;
  const double q8 = q5 * q5 + delt * q4 * q5 * t2;
  const double q9 = 1.0 + 2.0 * b * t2;
  const double hb = -bet * b * t6 * ( 2.0 + b * t2 ) / q8;
  const double hrs = -rsthrd * hb * bec * ecrs;
  const double ht = 2.0 * bet * q9 / q8;

  *ec += h;
  *vc1 = vc + h + hrs - t2 * ht * seven_sixth;
  *vc2 = -ht / ( rho * twoks * twoks );

}

// spin polarized case
void RSHFunctional::PBE_correlation_sp(const double rho_up, const double rho_dn,
  const double grad_up, const double grad_dn, const double grad, double *ec,
  double *vc1_up, double *vc1_dn, double *vc2)
{
  *ec = 0.0;
  *vc1_up = 0.0;
  *vc1_dn = 0.0;
  *vc2 = 0.0;

  const double rh_up = ( rho_up < 1.e-18 ) ? 0.0 : rho_up;
  const double rh_dn = ( rho_dn < 1.e-18 ) ? 0.0 : rho_dn;

  const double third = 1.0 / 3.0;
  const double third2 = 2.0 / 3.0;
  const double third4 = 4.0 / 3.0;
  const double sixthm = -1.0 / 6.0;
  const double pi32third = 3.09366772628014; /* (3*pi^2 ) ^(1/3) */
  const double alpha = 1.91915829267751; /* pow(9.0*pi/4.0, third)*/
  const double seven_sixth = 7.0 / 6.0;
  const double four_over_pi = 1.27323954473516;
  const double gam = 0.5198420997897463; /* gam = 2^(4/3) - 2 */
  const double fzz = 8.0 / ( 9.0 * gam );
  const double gamma = 0.03109069086965489; /* gamma = (1-ln2)/pi^2 */
  const double bet = 0.06672455060314922; /* see [a] (4) */
  const double delt = bet / gamma;
  const double eta = 1.e-12; // small number to avoid blowup as |zeta|->1

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

  const double rhotot = rh_up + rh_dn;

  const double rh13 = pow(rhotot,third);
  const double zet = ( rh_up - rh_dn ) / rhotot;
  const double g = 0.5 * ( pow(1.0 + zet,third2) + pow(1.0 - zet,third2) );
  const double fk = pi32third * rh13;
  const double rs = alpha / fk;
  const double twoksg = 2.0 * sqrt(four_over_pi * fk) * g;
  const double t = grad / ( twoksg * rhotot );

  const double rtrs = sqrt(rs);
  double eu, eurs, ep, eprs, alfm, alfrsm;
  gcor2(0.0310907,0.2137,7.5957,3.5876,1.6382,0.49294,rtrs,&eu,&eurs);
  gcor2(0.01554535,0.20548,14.1189,6.1977,3.3662,0.62517,rtrs,&ep,&eprs);
  gcor2(0.0168869,0.11125,10.357,3.6231,0.88026,0.49671,rtrs,&alfm,&alfrsm);
  const double z4 = zet * zet * zet * zet;
  const double f = ( pow(1.0 + zet,third4) + pow(1.0 - zet,third4) - 2.0 )
    / gam;
  *ec = eu * ( 1.0 - f * z4 ) + ep * f * z4 - alfm * f * ( 1.0 - z4 ) / fzz;

  /* LSD potential from [c] (A1) */
  /* ecrs = d ec / d rs [c] (A2) */
  const double ecrs = eurs * ( 1.0 - f * z4 ) + eprs * f * z4 - alfrsm * f
    * ( 1.0 - z4 ) / fzz;
  const double fz = third4 * ( pow(1.0 + zet,third) - pow(1.0 - zet,third) )
    / gam;
  const double eczet = 4.0 * ( zet * zet * zet ) * f * ( ep - eu + alfm / fzz )
    + fz * ( z4 * ep - z4 * eu - ( 1.0 - z4 ) * alfm / fzz );
  const double comm = *ec - rs * ecrs * third - zet * eczet;
  *vc1_up = comm + eczet;
  *vc1_dn = comm - eczet;

  /* PBE correlation energy */
  /* b = A of [a] (8) */

  const double g3 = g * g * g;
  const double pon = -*ec / ( g3 * gamma );
  const double b = delt / ( exp(pon) - 1.0 );
  const double b2 = b * b;
  const double t2 = t * t;
  const double t4 = t2 * t2;
  const double q4 = 1.0 + b * t2;
  const double q5 = q4 + b2 * t4;
  const double h = g3 * gamma * log(1.0 + delt * q4 * t2 / q5);

  /* Energy done, now the potential, using appendix E of [b] */

  const double g4 = g3 * g;
  const double t6 = t4 * t2;
  const double rsthrd = rs * third;
  const double gz = ( pow(( 1.0 + zet ) * ( 1.0 + zet ) + eta,sixthm) - pow(
    ( 1.0 - zet ) * ( 1.0 - zet ) + eta,sixthm) ) * third;
  const double fac = delt / b + 1.0;
  const double bg = -3.0 * b2 * *ec * fac / ( bet * g4 );
  const double bec = b2 * fac / ( bet * g3 );
  const double q8 = q5 * q5 + delt * q4 * q5 * t2;
  const double q9 = 1.0 + 2.0 * b * t2;
  const double hb = -bet * g3 * b * t6 * ( 2.0 + b * t2 ) / q8;
  const double hrs = -rsthrd * hb * bec * ecrs;
  const double hzed = 3.0 * gz * h / g + hb * ( bg * gz + bec * eczet );
  const double ht = 2.0 * bet * g3 * q9 / q8;

  double ccomm = h + hrs - t2 * ht * seven_sixth;
  const double pref = hzed - gz * t2 * ht / g;

  ccomm -= pref * zet;

  *ec += h;

  *vc1_up += ccomm + pref;
  *vc1_dn += ccomm - pref;

  *vc2 = -ht / ( rhotot * twoksg * twoksg );

}

// update exchange correlation energy and potential
void RSHFunctional::setxc(void)
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert(rho != 0);
    assert(grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0);
    assert(exc != 0);
    assert(vxc1 != 0);
    assert(vxc2 != 0);

#pragma omp parallel for
    for ( int i = 0; i < _np; i++ )
    {
      // evaluate gradient
      const double grad = sqrt(grad_rho[0][i] * grad_rho[0][i] + grad_rho[1][i]
        * grad_rho[1][i] + grad_rho[2][i] * grad_rho[2][i]);

      // calculate HSE exchange and PBE correlation
      double ex, vx1, vx2, ec, vc1, vc2;
      RSH_exchange(rho[i],grad,1 - x_coeff_,omega,&ex,&vx1,&vx2);
      PBE_correlation(rho[i],grad,&ec,&vc1,&vc2);

      // combine exchange and correlation energy
      exc[i] = ex + c_coeff_ * ec;
      vxc1[i] = vx1 + c_coeff_ * vc1;
      vxc2[i] = vx2 + c_coeff_ * vc2;
    }

  }
  else
  {
    assert(rho_up != 0);
    assert(rho_dn != 0);
    assert(grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0);
    assert(grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0);
    assert(exc_up != 0);
    assert(exc_dn != 0);
    assert(vxc1_up != 0);
    assert(vxc1_dn != 0);
    assert(vxc2_upup != 0);
    assert(vxc2_updn != 0);
    assert(vxc2_dnup != 0);
    assert(vxc2_dndn != 0);

#pragma omp parallel for
    for ( int i = 0; i < _np; i++ )
    {
      // evaluate gradient
      double grx_up = grad_rho_up[0][i];
      double gry_up = grad_rho_up[1][i];
      double grz_up = grad_rho_up[2][i];
      double grx_dn = grad_rho_dn[0][i];
      double gry_dn = grad_rho_dn[1][i];
      double grz_dn = grad_rho_dn[2][i];
      double grx = grx_up + grx_dn;
      double gry = gry_up + gry_dn;
      double grz = grz_up + grz_dn;
      double grad_up =
        sqrt(grx_up * grx_up + gry_up * gry_up + grz_up * grz_up);
      double grad_dn =
        sqrt(grx_dn * grx_dn + gry_dn * gry_dn + grz_dn * grz_dn);
      double grad = sqrt(grx * grx + gry * gry + grz * grz);

      // calculate HSE exchange and PBE correlation
      double ex_up, vx1_up, vx2_up, ex_dn, vx1_dn, vx2_dn;
      double ec, vc1_up, vc1_dn, vc2;
      RSH_exchange(2.0 * rho_up[i],2.0 * grad_up,1 - x_coeff_,omega,&ex_up,
        &vx1_up,&vx2_up);
      RSH_exchange(2.0 * rho_dn[i],2.0 * grad_dn,1 - x_coeff_,omega,&ex_dn,
        &vx1_dn,&vx2_dn);
      PBE_correlation_sp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&ec,&vc1_up,
        &vc1_dn,&vc2);

      // combine exchange and correlation energy
      exc_up[i] = ex_up + c_coeff_ * ec;
      exc_dn[i] = ex_dn + c_coeff_ * ec;
      vxc1_up[i] = vx1_up + c_coeff_ * vc1_up;
      vxc1_dn[i] = vx1_dn + c_coeff_ * vc1_dn;
      vxc2_upup[i] = 2 * vx2_up + c_coeff_ * vc2;
      vxc2_dndn[i] = 2 * vx2_dn + c_coeff_ * vc2;
      vxc2_updn[i] = c_coeff_ * vc2;
      vxc2_dnup[i] = c_coeff_ * vc2;
    }
  }
}
