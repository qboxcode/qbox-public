////////////////////////////////////////////////////////////////////////////////
//
// sinft.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: sinft.C,v 1.2 2007-08-13 21:23:51 fgygi Exp $

#include "sinft.h"
#include <math.h>
#include <assert.h>
#include "fftw.h"
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
