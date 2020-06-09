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
// test_fftw.cpp
//
////////////////////////////////////////////////////////////////////////////////

#ifdef USE_FFTW2

#include "Timer.h"

#include <iostream>
#include <complex>
#include <valarray>
using namespace std;
#include <cassert>

#include "fftw.h"

#ifdef IA32
#include "readTSC.h"
long long clk, clk_bwd, clk_fwd;
#endif

int main(int argc, char**argv)
{
  const int niter = 10;
  const int np = atoi(argv[1]);
  const int nvec = atoi(argv[2]);
  const int ldz = np + 1;

  fftw_plan fwplan, bwplan;

  // resize array zvec holding columns
  valarray<complex<double> > zvec(nvec * ldz);
  //cout << "zvec ptr: " << &zvec[0] << endl;

  // initialization of FFT libs

// #define FFTWMEASURE 1
#if FFTWMEASURE
  // FFTWMEASURE
  fwplan = fftw_create_plan(np,FFTW_FORWARD,FFTW_MEASURE|FFTW_IN_PLACE);
  bwplan = fftw_create_plan(np,FFTW_BACKWARD,FFTW_MEASURE|FFTW_IN_PLACE);
#else
  // FFTW_ESTIMATE
  fwplan = fftw_create_plan(np,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
  bwplan = fftw_create_plan(np,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif

  Timer t_fwd,t_bwd;

#ifdef IA32
  clk_bwd = 0;
  clk_fwd = 0;
#endif
  for ( int iter = 0; iter < niter; iter++ )
  {
  t_bwd.start();

   /*
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
  int ntrans = nvec;
  int inc1 = 1;
  int inc2 = ldz;
#ifdef IA32
  clk = readTSC();
#endif
  fftw(bwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
#ifdef IA32
  clk_bwd += readTSC() - clk;
#endif
  t_bwd.stop();
  t_fwd.start();
#ifdef IA32
  clk = readTSC();
#endif
  fftw(fwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
#ifdef IA32
  clk_fwd += readTSC() - clk;
#endif
  t_fwd.stop();
  }

  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);

  cout << " fwd: time per transform (in-place,generic)"
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_fwd.real()/(niter*nvec) << " microseconds"
#ifdef IA32
       << "  " << clk_fwd/(niter*nvec) << " cycles"
#endif
       << endl;

  cout << " bwd: time per transform (in-place,generic)"
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_bwd.real()/(niter*nvec) << " microseconds"
#ifdef IA32
       << "  " << clk_bwd/(niter*nvec) << " cycles"
#endif
       << endl;

#if 1
  // Use out-of-place, specific plan

  valarray<complex<double> > zvec_out(zvec.size());
  t_bwd.reset();
  t_fwd.reset();

  fwplan = fftw_create_plan_specific(np,
    FFTW_FORWARD,FFTW_ESTIMATE|FFTW_OUT_OF_PLACE,
    (FFTW_COMPLEX*)&zvec[0],1,(FFTW_COMPLEX*)&zvec_out[0],1);
  bwplan = fftw_create_plan_specific(np,
    FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_OUT_OF_PLACE,
    (FFTW_COMPLEX*)&zvec[0],1,(FFTW_COMPLEX*)&zvec_out[0],1);

#ifdef IA32
  clk_bwd = 0;
  clk_fwd = 0;
#endif
  for ( int iter = 0; iter < niter; iter++ )
  {

    int ntrans = nvec;
    int inc1 = 1;
    int inc2 = ldz;
    t_bwd.start();
#ifdef IA32
  clk = readTSC();
#endif
    fftw(bwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                       (FFTW_COMPLEX*)&zvec_out[0],inc1,inc2);
#ifdef IA32
  clk_bwd += readTSC() - clk;
#endif
    t_bwd.stop();

    t_fwd.start();
#ifdef IA32
  clk = readTSC();
#endif
    fftw(fwplan,ntrans,(FFTW_COMPLEX*)&zvec[0],inc1,inc2,
                       (FFTW_COMPLEX*)&zvec_out[0],inc1,inc2);
#ifdef IA32
  clk_fwd += readTSC() - clk;
#endif
    t_fwd.stop();

  }

  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);

  cout << " fwd: time per transform (out-of-place,specific)"
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_fwd.real()/(niter*nvec) << " microseconds"
#ifdef IA32
       << "  " << clk_fwd/(niter*nvec) << " cycles"
#endif
       << endl;

  cout << " bwd: time per transform (out-of-place,specific)"
#if FFTWMEASURE
       << "(fftw-measure)"
#endif
       << ": " << 1.e6*t_bwd.real()/(niter*nvec) << " microseconds"
#ifdef IA32
       << "  " << clk_bwd/(niter*nvec) << " cycles"
#endif
       << endl;
#endif
  return 0;
}
#else
int main(int argc, char**argv)
{}
#endif // USE_FFTW2
