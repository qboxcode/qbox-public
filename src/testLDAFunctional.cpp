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
// testLDAFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////

// Test the LDA functional by computing the xc energy of a gaussian
// of width 0.1 a.u. in a cube of side 1.0 a.u.
// With a cube of side 1.0 and 32x32x32 points,
// The xc energy must be -2.8105 a.u.
// dExc/da must be 0.911682

#include<iostream>
#include<vector>
#include "LDAFunctional.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
  if ( argc != 3 )
  {
    cout << " use: testLDAFunctional alat n" << endl;
    return 1;
  }
  double a = atof(argv[1]);
  double omega = a*a*a;
  int n = atoi(argv[2]);
  int n3 = n*n*n;
  vector<vector<double> > rh(1);
  rh[0].resize(n3);
  double excsum = 0.0, dxcsum = 0.0;

  double rc = 0.1 * a;
  double pi = 4.0 * atan(1.0);
  double fac = 1.0 / ( pow(pi,1.5) * rc*rc*rc );
  double sum = 0.0;

  for ( int i = 0; i < n; i++ )
  {
    double x = ( i * a ) / n - a/2;
    for ( int j = 0; j < n; j++ )
    {
      double y = ( j * a ) / n - a/2;
      for ( int k = 0; k < n; k++ )
      {
        double z = ( k * a ) / n - a/2;
        double r2 = x*x + y*y + z*z;
        int ii = i + n * ( j + n * k );
        rh[0][ii] = fac * exp( -r2 / (rc*rc) );
        sum += rh[0][ii];
      }
    }
  }
  sum = sum * omega / n3;
  // the density should be normalized
  cout << " Integrated density: " << sum << endl;

  LDAFunctional xcf(rh);
  xcf.setxc();

  for ( int i = 0; i < n3; i++ )
    excsum += xcf.rho[i] * xcf.exc[i];
  for ( int i = 0; i < n3; i++ )
    dxcsum += xcf.rho[i] * ( xcf.exc[i] - xcf.vxc1[i] );

  cout << " Total LDA xc energy: " << excsum * omega / n3 << endl;

  // Note: the energy variation is 3 * dExc/da * delta(a)
  cout << " dExc/da: " << dxcsum * omega / ( n3 * a ) << endl;

  return 0;
}
