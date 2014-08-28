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

#include "sinft.h"
#include <math.h>
#include <assert.h>
#if USE_FFTW2
#if USE_DFFTW
#include "dfftw.h"
#else
#include "fftw.h"
#endif
#endif
#if USE_FFTW3
#include "fftw3.h"
#endif
#include <vector>
#include <complex>
using namespace std;

void sinft(int n, double *y)
{
  vector<complex<double> > zin(2*n), zout(2*n);
  fftw_plan fwplan;
#if USE_FFTW2
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
#elif USE_FFTW3
  fwplan = fftw_plan_dft_1d(2*n,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0],
                            FFTW_FORWARD, FFTW_ESTIMATE);
#else
#error
#endif

  zin[0] = 0.0;
  for ( int i = 1; i < n; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = -t;
  }
#if USE_FFTW2
  fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
#elif USE_FFTW3
  fftw_execute(fwplan);
#else
#error
#endif
  for ( int i = 0; i < n; i++ )
  {
    y[i] = -0.5 * imag(zout[i]);
  }
  fftw_destroy_plan(fwplan);
}

void cosft1(int n, double *y)
{
  /* Note: the array y contains n+1 elements */
  vector<complex<double> > zin(2*n), zout(2*n);
  fftw_plan fwplan;
#if USE_FFTW2
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
#elif USE_FFTW3
  fwplan = fftw_plan_dft_1d(2*n,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0],
                            FFTW_FORWARD, FFTW_ESTIMATE);
#else
#error
#endif

  zin[0] = y[0];
  for ( int i = 1; i < n+1; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = t;
  }
#if USE_FFTW2
  fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
#elif USE_FFTW3
  fftw_execute(fwplan);
#else
#error
#endif
  y[0] = 0.5 * real(zout[0]);
  for ( int i = 1; i < n; i++ )
  {
    y[i] = 0.5 * real(zout[i]);
  }
  fftw_destroy_plan(fwplan);
}
