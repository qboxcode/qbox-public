////////////////////////////////////////////////////////////////////////////////
//
// test_fftw.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Timer.h"

#include <iostream>
#include <complex>
#include <valarray>
using namespace std;
#include <cassert>

#include "fftw.h"
extern "C" {
  void zdscal_(int *,double *,complex<double> *,int *);
}

int main(int argc, char**argv)
{
  const int niter = 10;
  const int np2_ = atoi(argv[1]);
  const int nvec_ = atoi(argv[2]);
  const int ldz = np2_ + 1;

  fftw_plan fwplan2, bwplan2;

  // resize array zvec holding columns
  valarray<complex<double> > zvec_(nvec_ * ldz);
  
  // initialization of FFT libs

#if FFTWMEASURE
  // FFTWMEASURE
  fwplan2 = fftw_create_plan(np2_,FFTW_FORWARD,FFTW_MEASURE|FFTW_IN_PLACE);
  bwplan2 = fftw_create_plan(np2_,FFTW_BACKWARD,FFTW_MEASURE|FFTW_IN_PLACE);
#else
  // FFTW_ESTIMATE
  fwplan2 = fftw_create_plan(np2_,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
  bwplan2 = fftw_create_plan(np2_,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_IN_PLACE);
#endif

  Timer t_fwd,t_bwd;

  for ( int iter = 0; iter < niter; iter++ )
  {
  t_bwd.start();

   /* 
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
  int ntrans = nvec_;
  int inc1 = 1;
  int inc2 = ldz;
  fftw(bwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  t_bwd.stop();
  

  // transform along z
  
 /*
  * void fftw(fftw_plan plan, int howmany,
  *    FFTW_COMPLEX *in, int istride, int idist,
  *    FFTW_COMPLEX *out, int ostride, int odist);
  */
  t_fwd.start();
  fftw(fwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  t_fwd.stop();

  int len = zvec_.size();
  double fac = 3.14;
  zdscal_(&len,&fac,&zvec_[0],&inc1);
  
  }

  fftw_destroy_plan(fwplan2);
  fftw_destroy_plan(bwplan2);

  cout << " fwd: " << t_fwd.real()/niter << endl;
  cout << " fwd: time per transform (microseconds): " 
       << 1.e6*t_fwd.real()/(niter*nvec_)
       << endl;
  cout << " bwd: " << t_bwd.real()/niter << endl;
  cout << " bwd: time per transform (microseconds): " 
       << 1.e6*t_bwd.real()/(niter*nvec_)
       << endl;

  return 0;
}
