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
// sinft.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: sinft.C,v 1.4 2008-09-08 15:56:20 fgygi Exp $

#include "sinft.h"
#include <math.h>
#include <assert.h>
#if USE_FFTW
#if USE_DFFTW
#include "dfftw.h"
#else
#include "fftw.h"
#endif
#endif
#include <vector>
#include <complex>
using namespace std;

void sinft(int n, double *y)
{
  fftw_plan fwplan;
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
  vector<complex<double> > zin(2*n), zout(2*n);
  zin[0] = 0.0;
  for ( int i = 1; i < n; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = -t;
  }
  fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
  for ( int i = 0; i < n; i++ )
  {
    y[i] = -0.5 * imag(zout[i]);
  }
  fftw_destroy_plan(fwplan);
}

void cosft1(int n, double *y)
{
  /* Note: the array y contains n+1 elements */
  fftw_plan fwplan;
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
  vector<complex<double> > zin(2*n), zout(2*n);

  zin[0] = y[0];
  for ( int i = 1; i < n+1; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = t;
  }
  fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
  y[0] = 0.5 * real(zout[0]);
  for ( int i = 1; i < n; i++ )
  {
    y[i] = 0.5 * real(zout[i]);
  }
  fftw_destroy_plan(fwplan);
}
