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
// ExponentialIntegral.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Calculate exponential integral using the algorithm of
// Tseng, Lee, Journal of Hydrology, 205 (1998) 38-51
//
////////////////////////////////////////////////////////////////////////////////

#include "ExponentialIntegral.h"
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

namespace util
{

// Calculate the exponential integral E_1(x):
//
//        inf
//          /     -t
//         |     e
// E (x) = | dt -----
//  1      |      t
//        /
//         x
//
// Input:  x - position at which exponential integral is evaluated (x > 0)
// Return: E1(x)
double E1(const double x)
{

  if (x < series_cutoff)
  {
    // use series expansion
    return E1_series(x);
  }
  else
  {
    // use gauss_laguerre expression
    return exp(-x) * gauss_laguerre(x);
  }
}

// Series expansion of the exponential integral
//
//                          n_cut
//                          -----     n  n
//                           \    (-1)  x
// E (x) = -gamma - ln(x) -   )   --------
//  1                        /     n * n!
//                          -----
//                          n = 1
//
// where gamma is the Euler constant.
// n_cut is set to 25
// Input:  x - position at which exponential integral is evaluated (x > 0)
// Return: approximation by series expansion for E_1(x)
double E1_series(const double x)
{
  // Euler constant
  const double EULER_GAMMA = 0.57721566490153286060651209008241;

  // Cutoff for series expansion
  const int itmax = 25;

  // initialize summation result
  double res = 0.0;

  // perform the summation
  for (int it = itmax; it > 1; it--)
  {
    // calculate 1/n
    const double fact = 1.0 / it;
    // add next term of summation
    res = x * fact * (fact - res);
  }

  // add everything up
  return -EULER_GAMMA - log(x) + x * (1.0 - res);
}

// The Gauss Laguerre expansion of the exponential integral can be written as
//
//             N
// E (x0)    -----     a
//  1         \         n
// ------ =    )   ---------
//   -x0      /     x  + x0
//  e        -----   n
//            n=1
//
// where the a_n and x_n are determined by least quadrature (see reference)
// Input: x0 - point at which Gaussian Laguerre quadrature is calculated
// Return: E_1(x0) / exp(-x0) in this approximation
double gauss_laguerre(const double x0)
{
  // initialize constants a_n and x_n
  const double size = 15;
  const double a[] = { 0.2182348859400869e+00, 0.3422101779228833e+00,
      0.2630275779416801e+00, 0.1264258181059305e+00, 0.4020686492100091e-01,
      0.8563877803611838e-02, 0.1212436147214252e-02, 0.1116743923442519e-03,
      0.6459926762022901e-05, 0.2226316907096273e-06, 0.4227430384979365e-08,
      0.3921897267041089e-10, 0.1456515264073126e-12, 0.1483027051113301e-15,
      0.1600594906211133e-19 };
  const double x[] = { 0.9330781201728180e-01, 0.4926917403018839e+00,
      0.1215595412070949e+01, 0.2269949526203743e+01, 0.3667622721751437e+01,
      0.5425336627413553e+01, 0.7565916226613068e+01, 0.1012022856801911e+02,
      0.1313028248217572e+02, 0.1665440770832996e+02, 0.2077647889944877e+02,
      0.2562389422672878e+02, 0.3140751916975394e+02, 0.3853068330648601e+02,
      0.4802608557268579e+02 };

  // initialize
  double res = 0.0;

  // evaluate a_n / ( x_n + x0 )
  for ( int i = 0; i < size; i++ )
  {
    res += a[i] / (x[i] + x0);
  }

  return res;
}

}
