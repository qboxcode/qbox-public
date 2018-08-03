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
// testFunction3d.C
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include "Function3d.h"
using namespace std;

// usage: ./testFunction3d file.xml

int main(int argc, char** argv)
{
  Function3d f;
  if ( argc != 2 )
  {
    cerr << "usage: " << argv[0] << " file.xml" << endl;
    return 1;
  }
  f.read(argv[1]);

  cout << "function name: " << f.name << endl;
  cout << "function grid size: " << f.nx << " " << f.ny << " " << f.nz  << endl;
  cout << "function array size: " << f.val.size() << endl;

  double sum = 0.0;
  for ( int i = 0; i < f.val.size(); i++ )
  {
    sum += f.val[i]*f.val[i];
  }
  cout << "function norm2: " << sum / f.val.size() << endl;

  // modify function
  f.name = "output";
  for ( int i = 0; i < f.val.size(); i++ )
  {
    f.val[i] *= 0.5;
  }

  f.write("output.xml");
}
