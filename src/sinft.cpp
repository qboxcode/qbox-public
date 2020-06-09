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
// sinft.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "sinft.h"
#include <math.h>
#include <assert.h>
#include <vector>
#include <complex>
using namespace std;

#if USE_FFTW2
#if USE_DFFTW
#include "dfftw.h"
#else
#include "fftw.h"
#endif
#elif USE_FFTW3
#include "fftw3.h"
#elif USE_ESSL_FFT
extern "C" {
  void dcft_(int *initflag, std::complex<double> *x, int *inc2x, int *inc3x,
             std::complex<double> *y, int *inc2y, int *inc3y,
             int *length, int *ntrans, int *isign,
             double *scale, double *aux1, int *naux1,
             double *aux2, int *naux2);
}
#else
// no FFT library
void cfft ( int idir, complex<double> *z1, complex<double> *z2, int n,
  int *inzee );
void fftstp ( int idir, complex<double> *zin, int after,
              int now, int before, complex<double> *zout );
#endif

void sinft(int n, double *y)
{
  vector<complex<double> > zin(2*n), zout(2*n);
#if defined(USE_FFTW2) || defined(USE_FFTW3)
  fftw_plan fwplan;
#endif
#if USE_FFTW2
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
#elif USE_FFTW3
  fwplan = fftw_plan_dft_1d(2*n,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0],
                            FFTW_FORWARD, FFTW_ESTIMATE);
#elif USE_ESSL_FFT
  int np = 2 * n;
  int naux1 = (int) (30000 + 2.28 * np);
  std::vector<double> aux1(naux1);
  int ntrans = 1;
  int naux2 = (int) (20000 + 2.28 * np + (256 + 2*np)*min(64,ntrans));
  std::vector<double> aux2(naux2);
#else
  // no FFT library
  // no initialization needed
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
#elif USE_ESSL_FFT
  // initialize forward transform
  int initflag = 1;
  int inc1 = 1, inc2 = np;
  int isign = 1;
  double scale = 1.0;
  complex<double> *p = 0;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np,&ntrans,
        &isign,&scale,&aux1[0],&naux1,&aux2[0],&naux2);

  // call transform
  initflag = 0;
  dcft_(&initflag,&zin[0],&inc1,&inc2,&zout[0],&inc1,&inc2,&np,&ntrans,
        &isign,&scale,&aux1[0],&naux1,&aux2[0],&naux2);
#else
  // no FFT library
  int idir = 1, inzee = 1, np = 2*n;
  cfft ( idir, &zin[0],&zout[0],np, &inzee );
#endif
  for ( int i = 0; i < n; i++ )
  {
    y[i] = -0.5 * imag(zout[i]);
  }
#if defined(USE_FFTW2) || defined(USE_FFTW3)
  fftw_destroy_plan(fwplan);
#endif
}

void cosft1(int n, double *y)
{
  /* Note: the array y contains n+1 elements */
  vector<complex<double> > zin(2*n), zout(2*n);
#if defined(USE_FFTW2) || defined(USE_FFTW3)
  fftw_plan fwplan;
#endif
#if USE_FFTW2
  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
#elif USE_FFTW3
  fwplan = fftw_plan_dft_1d(2*n,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0],
                            FFTW_FORWARD, FFTW_ESTIMATE);
#elif USE_ESSL_FFT
  int np = 2 * n;
  int naux1 = (int) (30000 + 2.28 * np);
  std::vector<double> aux1(naux1);
  int ntrans = 1;
  int naux2 = (int) (20000 + 2.28 * np + (256 + 2*np)*min(64,ntrans));
  std::vector<double> aux2(naux2);
#else
  // no FFT library
  // no initialization needed
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
#elif USE_ESSL_FFT
  // initialize forward transform
  int initflag = 1;
  int inc1 = 1, inc2 = np;
  int isign = 1;
  double scale = 1.0;
  complex<double> *p = 0;
  dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np,&ntrans,
        &isign,&scale,&aux1[0],&naux1,&aux2[0],&naux2);

  // call transform
  initflag = 0;
  dcft_(&initflag,&zin[0],&inc1,&inc2,&zout[0],&inc1,&inc2,&np,&ntrans,
        &isign,&scale,&aux1[0],&naux1,&aux2[0],&naux2);
#else
  // no FFT library
  int idir = 1, inzee = 1, np = 2*n;
  cfft ( idir, &zin[0],&zout[0],np, &inzee );
#endif
  y[0] = 0.5 * real(zout[0]);
  for ( int i = 1; i < n; i++ )
  {
    y[i] = 0.5 * real(zout[i]);
  }
#if defined(USE_FFTW2) || defined(USE_FFTW3)
  fftw_destroy_plan(fwplan);
#endif
}
