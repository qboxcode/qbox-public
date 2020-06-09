////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 The Regents of the University of California
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
// test_vext.cpp
//
// Generate an external potential in XML format for input to the response
// command.
// vext = sine wave in the x direction
// use: test_vext a b c np0 np1 np2
//
////////////////////////////////////////////////////////////////////////////////

#include "Function3d.h"
#include<cassert>
#include<vector>
#include<string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  if ( argc == 1 )
  {
    cerr << "use: test_vext a b c np0 np1 np2" << endl;
    return 1;
  }
  double a = atof(argv[1]);
  double b = atof(argv[2]);
  double c = atof(argv[3]);
  Function3d f;
  f.a = D3vector(a,0,0);
  f.b = D3vector(0,b,0);
  f.c = D3vector(0,0,c);
  f.nx = atoi(argv[4]);
  f.ny = atoi(argv[5]);
  f.nz = atoi(argv[6]);
  f.val.resize(f.nx*f.ny*f.nz);
  f.name = "delta_v";
  for ( int i = 0; i < f.nx; i++ )
    for ( int j = 0; j < f.ny; j++ )
      for ( int k = 0; k < f.nz; k++ )
      {
        double x = ( a * i ) / f.nx;
        double y = ( b * j ) / f.ny;
        double z = ( c * k ) / f.nz;
        f.val[i+f.nx*(j+f.ny*k)] = sin(2*M_PI*x/a);
      }
  cout << f;
  return 0;
}
