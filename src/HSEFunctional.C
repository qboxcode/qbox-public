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
#include <cassert>
#include <cmath>

using namespace std;

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
void HSE_enhance(double s, double kF, double w, double *fx, double *dfx_ds,
    double* dfx_dkf) {

  // definition of constants
  const double A = 1.0161144, B = -0.37170836, C = -0.077215461,
      D = 0.57786348, E = -0.051955731;
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

}


void HSEFunctional::setxc(void) {
  // dummy
}
