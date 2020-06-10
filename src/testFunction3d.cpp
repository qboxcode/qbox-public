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
// testFunction3d.cpp
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<iomanip>
#include "Function3d.h"
using namespace std;

// usage: ./testFunction3d file.xml

// Read a Function3d from an XML file and print contents as a cube file

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

  f.write("output.xml");

  // output data in cube format on stdout
  cout << "Created by testFunction3d" << endl;
  cout << endl;
  cout << "8 0 0 0" << endl;
  cout << f.nx << " " << f.a / f.nx << endl;
  cout << f.ny << " " << f.b / f.ny << endl;
  cout << f.nz << " " << f.c / f.nz << endl;
  cout << "1 1 0.0 0.0 0.0" << endl;
  cout << "1 1 " << f.a << endl;
  cout << "1 1 " << f.b << endl;
  cout << "1 1 " << f.c << endl;
  cout << "1 1 " << f.a+f.b << endl;
  cout << "1 1 " << f.b+f.c << endl;
  cout << "1 1 " << f.c+f.a << endl;
  cout << "1 1 " << f.a+f.b+f.c << endl;
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(5);
  for ( int i = 0; i < f.val.size(); i++ )
  {
    cout << setw(13) << f.val[i];
    if ( ( (i+1) % 6 ) == 0 )
      cout << endl;
  }
  cout << endl;
}
