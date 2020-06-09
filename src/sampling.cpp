////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2010 The Regents of the University of California
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
// sampling.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <cmath> // sqrt, log
#include <cstdlib> // drand48
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void normal_dev(double* d1, double *d2)
{
  // Box-Muller algorithm: generate a pair of normal random deviates
  // J. R. Bell, Algorithm 334, Comm. ACM 11, 498, (1968).
  // modified following R. Knop, Comm. ACM 12, 281 (1969).
  double x,y,s;
  do
  {
    x = 2.0 * drand48() - 1.0;
    y = 2.0 * drand48() - 1.0;
    s = x*x + y*y;
  }
  while ( s > 1.0 );
  double l = sqrt(-2.0 * log(s)/s);
  *d1 = x * l;
  *d2 = y * l;
}

////////////////////////////////////////////////////////////////////////////////
// generate deviates for the gamma function
// use the sum of nsum gaussian deviates
// use: tgamma_sum nsum nval seed
////////////////////////////////////////////////////////////////////////////////
double gamma_dev(int a)
{
  double am=a-1.0, e, s=sqrt(2.0*(a-1.0)+1.0), v1, v2, x, y;
  do
  {
    do
    {
      do
      {
        v1 = 2.0*drand48()-1.0;
        v2 = 2.0*drand48()-1.0;
      } while ( v1*v1 + v2*v2 > 1.0 || v1 == 0.0 );
      y = v2/v1;
      x = s*y+am;
    } while ( x <= 0.0 );
    e = (1.0+y*y)*exp(am*log(x/am)-s*y);
  } while ( drand48() > e );
  return x;
}
