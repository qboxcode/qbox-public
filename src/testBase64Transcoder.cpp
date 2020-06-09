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
// testBase64Transcoder.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Base64Transcoder.h"
#include <iostream>
#include <cassert>
using namespace std;

int main()
{
  const int n = 7;
  double a[n] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  double c[n] = { 0, 0, 0, 0, 0, 0, 0 };

  Base64Transcoder xcdr;

  int nbytes = n * sizeof(double);

  int nchars = xcdr.nchars(nbytes);

  cout << " nbytes=" << nbytes << endl;
  cout << " nchars=" << nchars << endl;

  char* b = new char[nchars];

  xcdr.encode(nbytes,(unsigned char*) &a[0],b);

  cout << " b=" << b << endl;

  xcdr.decode(nchars,b,(unsigned char*) &c[0]);

  for ( int i = 0; i < n; i++ )
    assert(a[i]==c[i]);

  cout << " done" << endl;

  return 0;
}
