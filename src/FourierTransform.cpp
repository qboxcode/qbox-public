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
// FourierTransform.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "FourierTransform.h"
#include "Basis.h"
#include "BasisMapping.h"
#include "blas.h"

#include <complex>
#include <cassert>
#include <cstdlib> // abs

#if _OPENMP
#include <omp.h>
#else
// _OPENMP is not defined
#if defined(USE_FFTW3_THREADS)
#error "Need OpenMP to use FFTW3 threads"
#endif
#endif

#if defined(USE_FFTW2) || defined(USE_FFTW3)
#ifdef ADD_
#define zdscal zdscal_
#define zcopy zcopy_
#endif
#endif

#if defined(USE_FFTW2)
#if defined(FFTWMEASURE)
#define FFTW_ALGO FFTW_MEASURE
#else
#define FFTW_ALGO FFTW_ESTIMATE
#endif
#endif

#if defined(USE_FFTW3)
#if defined(FFTWMEASURE)
#define FFTW_ALGO ( FFTW_MEASURE | FFTW_UNALIGNED )
#else
#define FFTW_ALGO ( FFTW_ESTIMATE | FFTW_UNALIGNED )
#endif
#endif

#if defined(USE_FFTW2) || defined(USE_FFTW3)
extern "C" void zdscal(int *n,double *alpha,std::complex<double> *x,int *incx);
#elif USE_ESSL_FFT
extern "C" {
  void dcft_(int *initflag, std::complex<double> *x, int *inc2x, int *inc3x,
             std::complex<double> *y, int *inc2y, int *inc3y,
             int *length, int *ntrans, int *isign,
             double *scale, double *aux1, int *naux1,
             double *aux2, int *naux2);
  void dcft2_(int *initflag, std::complex<double> *x, int *inc1x, int *inc2x,
             std::complex<double> *y, int *inc1y, int *inc2y,
             int *n1, int *n2, int *isign,
             double *scale, double *aux1, int *naux1,
             double *aux2, int *naux2);
#define USE_GATHER_SCATTER 1
}
#elif defined(FFT_NOLIB)
void cfftm ( std::complex<double> *ain, std::complex<double> *aout,
  double scale, int ntrans, int length, int ainc, int ajmp, int idir );
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if USE_GATHER_SCATTER
extern "C" {
  // zgthr: x(i) = y(indx(i))
  void zgthr_(int* n, std::complex<double>* y,
              std::complex<double>* x, int*indx);
  // zsctr: y(indx(i)) = x(i)
  void zsctr_(int* n, std::complex<double>* x, int* indx,
              std::complex<double>* y);
}
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////////////
FourierTransform::FourierTransform (const Basis &basis,
  int np0, int np1, int np2) : basis_(basis),
  bm_(BasisMapping(basis,np0,np1,np2)), np0_(np0), np1_(np1), np2_(np2),
  nvec_(bm_.nvec())
{
  // compute number of transforms along the x,y and z directions
  // ntrans0_ is the number of transforms along x in one of the two blocks
  // of vectors corresponding to positive and negative y indices
  ntrans0_ = max(abs(basis_.idxmax(1)),abs(basis_.idxmin(1)))+1;
  ntrans1_ = np0_;
  ntrans2_ = nvec_;

  // resize array zvec holding columns
  zvec_.resize(nvec_ * np2_);

  // Initialize FT library auxiliary arrays
  init_lib();
}

////////////////////////////////////////////////////////////////////////////////
FourierTransform::~FourierTransform()
{
#if USE_FFTW2
  fftw_destroy_plan(fwplan0);
  fftw_destroy_plan(fwplan1);
  fftw_destroy_plan(fwplan2);
  fftw_destroy_plan(bwplan0);
  fftw_destroy_plan(bwplan1);
  fftw_destroy_plan(bwplan2);
#endif

#if USE_FFTW3
#if USE_FFTW3_THREADS
  fftw_cleanup_threads();
#endif
#if defined(USE_FFTW3_2D) || defined(USE_FFTW3_THREADS)
  fftw_destroy_plan(fwplan2d);
  fftw_destroy_plan(bwplan2d);
#else
  fftw_destroy_plan(fwplanx);
  fftw_destroy_plan(bwplanx);
  fftw_destroy_plan(fwplany);
  fftw_destroy_plan(bwplany);
#endif
  fftw_destroy_plan(fwplan);
  fftw_destroy_plan(bwplan);
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::forward(complex<double>* f, complex<double>* c)
{
  fwd(f);
#if TIMING
  tm_map_fwd.start();
#endif
  bm_.zvec_to_vector(&zvec_[0],c);
#if TIMING
  tm_map_fwd.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::backward(const complex<double>* c, complex<double>* f)
{
#if TIMING
  tm_map_bwd.start();
#endif
  bm_.vector_to_zvec(c,&zvec_[0]);
#if TIMING
  tm_map_bwd.stop();
#endif
  bwd(f);
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::forward(complex<double>* f,
  complex<double>* c1, complex<double>* c2)
{
  fwd(f);
#if TIMING
  tm_map_fwd.start();
#endif
  bm_.zvec_to_doublevector(&zvec_[0],c1,c2);
#if TIMING
  tm_map_fwd.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::backward(const complex<double>* c1,
                               const complex<double>* c2,
                               complex<double>* f)
{
#if TIMING
  tm_map_bwd.start();
#endif
  bm_.doublevector_to_zvec(c1,c2,&zvec_[0]);
#if TIMING
  tm_map_bwd.stop();
#endif
  bwd(f);
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fwd(complex<double>* val)
{
#if TIMING
  tm_fwd.start();
#endif
  fxy(val);
#if TIMING
  tm_trans_fwd.start();
#endif
  bm_.transpose_fwd(val,&zvec_[0]);
#if TIMING
  tm_trans_fwd.stop();
#endif
  fz();
#if TIMING
  tm_fwd.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::bwd(complex<double>* val)
{
#if TIMING
  tm_bwd.start();
#endif
  fz_inv();
#if TIMING
  tm_trans_bwd.start();
#endif
  bm_.transpose_bwd(&zvec_[0],val);
#if TIMING
  tm_trans_bwd.stop();
#endif
  fxy_inv(val);
#if TIMING
  tm_bwd.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fxy(complex<double>* val)
{
#if TIMING
  tm_fxy.start();
#endif

//fftw_execute_dft is thread safe
#if USE_FFTW3
#if USE_FFTW3_THREADS
  fftw_execute_dft ( fwplan2d, (fftw_complex*)&val[0],
                     (fftw_complex*)&val[0] );
#elif USE_FFTW3_2D // USE_FFTW3_2D
  #pragma omp parallel for
  for ( int k = 0; k < np2_loc(); k++ )
    fftw_execute_dft ( fwplan2d, (fftw_complex*)&val[k*np0_*np1_],
                       (fftw_complex*)&val[k*np0_*np1_] );
#else // USE_FFTW3_2D
  for ( int k = 0; k < np2_loc(); k++ )
  {
    const int ibase = k * np0_ * np1_;
#if FFTW_TRANSPOSE
    #pragma omp parallel
    {
      vector<complex<double> >t_trans(np1_);
      #pragma omp for
      for ( int i = 0; i < np0_; i++ )
      {
        int length = t_trans.size();
        int inc1 = 1, inc2 = np0_;
        zcopy(&length, &val[ibase+i], &inc2, &t_trans[0], &inc1);
        fftw_execute_dft ( fwplany, (fftw_complex*)&t_trans[0],
                         (fftw_complex*)&t_trans[0]);
        zcopy(&length, &t_trans[0], &inc1, &val[ibase+i], &inc2);
      }
    }
#else // FFTW_TRANSPOSE
    #pragma omp parallel for
    for ( int i = 0; i < np0_; i++ )
    {
      fftw_execute_dft ( fwplany, (fftw_complex*)&val[ibase+i],
                         (fftw_complex*)&val[ibase+i]);
    }
#endif // FFTW_TRANSPOSE
    #pragma omp parallel for
    for ( int i = 0; i < ntrans0_; i++ )
    {
      // Transform first block along x: positive y indices
      fftw_execute_dft ( fwplanx,(fftw_complex*)&val[ibase+i*np0_],
                         (fftw_complex*)&val[ibase+i*np0_]);

      // Transform second block along x: negative y indices
      fftw_execute_dft ( fwplanx,
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_]);
    }
  }
#endif // USE_FFTW3_2D
#elif USE_ESSL_FFT
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if USE_ESSL_2DFFT

    // use 2D FFT for x and y transforms
    int inc1, inc2, istart, isign = 1, initflag = 0;
    double scale = 1.0;

    // xy transform
    istart = k * np0_ * np1_;
    inc1 = 1; inc2 = np0_;
    dcft2_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
          &np0_,&np1_,&isign,&scale,&aux1xyf[0],&naux1xy,&aux2[0],&naux2);

#else

    // use multiple 1-d FFTs for x and y transforms

    int inc1, inc2, ntrans, istart, length, isign = 1, initflag = 0;
    double scale = 1.0;
    // transform along y for all values of x
    ntrans = np0_;
    if ( ntrans > 0 )
    {
      inc1 = np0_;
      inc2 = 1;
      istart = k * np0_ * np1_;
      length = np1_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1yf[0],&naux1y,&aux2[0],&naux2);
    }

    // transform only non-zero vectors along x
    ntrans = ntrans0_;
    if ( ntrans > 0 )
    {
      inc1 = 1;
      inc2 = np0_;
      istart = k * np0_ * np1_;
      length = np0_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);

      inc1 = 1;
      inc2 = np0_;
      istart = np0_ * ( (np1_-ntrans) + k * np1_ );
      length = np0_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);
    }
#endif // USE_ESSL_2DFFT
  } // k
#elif USE_FFTW2
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if _OPENMP
  int ibase = k * np0_ * np1_;
  //complex<double> *tmp1 = new complex<double>[np1_];
  #pragma omp parallel for
  for ( int i = 0; i < np0_; i++ )
  {
    //#pragma omp task
    {
      // transform along y for all values of x
      // copy data to local array
      int one=1;
      #if 0
      zcopy_(&np1_,&val[ibase+i],&np0_,tmp1,&one);
      fftw_one(fwplan1,(FFTW_COMPLEX*)tmp1,(FFTW_COMPLEX*)0);
      zcopy_(&np1_,tmp1,&one,&val[ibase+i],&np0_);
      #else
      fftw(fwplan1,1,(FFTW_COMPLEX*)&val[ibase+i],np0_,one,
                     (FFTW_COMPLEX*)0,0,0);
      #endif
    }
  }
  //delete [] tmp1;

  #pragma omp parallel for
  for ( int i = 0; i < ntrans0_; i++ )
  {
    //#pragma omp task
    {
      // Transform first block along x: positive y indices
      fftw_one(fwplan0,(FFTW_COMPLEX*)&val[ibase+i*np0_],(FFTW_COMPLEX*)0);
      // Transform second block along x: negative y indices
      fftw_one(fwplan0,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                       (FFTW_COMPLEX*)0);
    }
  }
#else // _OPENMP
    int inc1, inc2, istart;

    // transform along y for all values of x
    int ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    fftw(fwplan1,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);

    ntrans = ntrans0_;
    // Transform first block along x: positive y indices
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    fftw(fwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
    // Transform second block along x: negative y indices
    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    fftw(fwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP
  } // k
#elif defined(FFT_NOLIB)
  // No library
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
    // transform along y for all values of x
    int ntrans = np0_;
    int istart = k * np0_ * np1_;
    int length = np1_;
    int ainc = np0_;
    int ajmp = 1;
    double scale = 1.0;
    int idir = 1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // transform along x for non-zero elements
    ntrans = ntrans0_;
    istart = k * np0_ * np1_;
    length = np0_;
    ainc   = 1;
    ajmp   = np0_;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );
  } // for k
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_fxy.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fxy_inv(complex<double>* val)
{
#if TIMING
  tm_fxy_inv.start();
#endif
#if USE_FFTW3
#if USE_FFTW3_THREADS
  fftw_execute_dft ( bwplan2d, (fftw_complex*)&val[0],
                     (fftw_complex*)&val[0] );
#elif USE_FFTW3_2D
  #pragma omp parallel for
  for ( int k = 0; k < np2_loc(); k++ )
    fftw_execute_dft ( bwplan2d, (fftw_complex*)&val[k*np0_*np1_],
                       (fftw_complex*)&val[k*np0_*np1_] );
#else // FFTW3_2D
  // fftw3 1d
  for ( int k = 0; k < np2_loc(); k++ )
  {
    int ibase = k * np0_ * np1_;
    #pragma omp parallel for
    for ( int i = 0; i < ntrans0_; i++ )
    {
      // Transform first block along x: positive y indices
      fftw_execute_dft ( bwplanx, (fftw_complex*)&val[ibase+i*np0_],
                         (fftw_complex*)&val[ibase+i*np0_]);
      // Transform second block along x: negative y indices
      fftw_execute_dft ( bwplanx,
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                         (fftw_complex*)&val[ibase+(np1_-ntrans0_+i)*np0_]);
    }
#if FFTW_TRANSPOSE
    #pragma omp parallel
    {
      vector<complex<double> >t_trans(np1_);
      #pragma omp for
      for ( int i = 0; i < np0_; i++ )
      {
        int length = t_trans.size();
        int inc1 = 1, inc2 = np0_;
        zcopy(&length, &val[ibase+i], &inc2, &t_trans[0], &inc1);
        fftw_execute_dft ( bwplany, (fftw_complex*)&t_trans[0],
                           (fftw_complex*)&t_trans[0]);
        zcopy(&length, &t_trans[0], &inc1, &val[ibase+i], &inc2);
      }
    }
#else // FFTW_TRANSPOSE
    #pragma omp parallel for
    for ( int i = 0; i < np0_; i++ )
    {
      fftw_execute_dft ( bwplany, (fftw_complex*)&val[ibase+i],
                         (fftw_complex*)&val[ibase+i]);
    }
#endif // FFTW_TRANSPOSE
  }
#endif // USE_FFTW3_2D

#elif USE_ESSL_FFT
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if USE_ESSL_2DFFT

    // use 2D FFT for x and y transforms
    int inc1, inc2, istart, isign = -1, initflag = 0;
    double scale = 1.0;

    // xy transform
    istart = k * np0_ * np1_;
    inc1 = 1; inc2 = np0_;
    dcft2_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
          &np0_,&np1_,&isign,&scale,&aux1xyb[0],&naux1xy,&aux2[0],&naux2);

#else

    // use multiple 1-d FFTs for x and y transforms

    int inc1, inc2, ntrans, istart, length, isign = -1, initflag = 0;
    double scale = 1.0;
    // transform only non-zero vectors along x
    // First block: positive y indices: [0,ntrans0_]
    ntrans = ntrans0_;
    if ( ntrans > 0 )
    {
      inc1 = 1;
      inc2 = np0_;
      istart = k * np0_ * np1_;
      length = np0_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);

      // Second block: negative y indices: [np1-ntrans0_,np1-1]
      inc1 = 1;
      inc2 = np0_;
      istart = np0_ * ( (np1_-ntrans) + k * np1_ );
      length = np0_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);
    }

    // transform along y for all values of x
    ntrans = np0_;
    if ( ntrans > 0 )
    {
      inc1 = np0_;
      inc2 = 1;
      istart = k * np0_ * np1_;
      length = np1_;
      dcft_(&initflag,&val[istart],&inc1,&inc2,&val[istart],&inc1,&inc2,
            &length,&ntrans,&isign,&scale,&aux1yb[0],&naux1y,&aux2[0],&naux2);
    }
#endif // USE_ESSL_2DFFT
  } // k

#elif USE_FFTW2
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
#if _OPENMP
  int ibase = k * np0_ * np1_;
  #pragma omp parallel for
  for ( int i = 0; i < ntrans0_; i++ )
  {
    //#pragma omp task
    {
      // Transform first block along x: positive y indices
      fftw_one(bwplan0,(FFTW_COMPLEX*)&val[ibase+i*np0_],(FFTW_COMPLEX*)0);
      //fftw(bwplan0,1,(FFTW_COMPLEX*)&val[ibase+i*np0_],1,np0_,
      //               (FFTW_COMPLEX*)0,0,0);
      // Transform second block along x: negative y indices
      fftw_one(bwplan0,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],
                       (FFTW_COMPLEX*)0);
      //fftw(bwplan0,1,(FFTW_COMPLEX*)&val[ibase+(np1_-ntrans0_+i)*np0_],1,np0_,
      //               (FFTW_COMPLEX*)0,0,0);
    }
  }

  //complex<double> *tmp1 = new complex<double>[np1_];
  #pragma omp parallel for
  for ( int i = 0; i < np0_; i++ )
  {
    {
      // transform along y for all values of x
      // copy data to local array
      int one=1;
      #if 0
      zcopy_(&np1_,&val[ibase+i],&np0_,tmp1,&one);
      fftw_one(bwplan1,(FFTW_COMPLEX*)tmp1,(FFTW_COMPLEX*)0);
      zcopy_(&np1_,tmp1,&one,&val[ibase+i],&np0_);
      #else
      fftw(bwplan1,1,(FFTW_COMPLEX*)&val[ibase+i],np0_,one,
                     (FFTW_COMPLEX*)0,0,0);
      #endif
    }
  }
  //delete [] tmp1;
#else // _OPENMP
    int inc1, inc2, istart;

    int ntrans = ntrans0_;
    // Transform first block along x: positive y indices
    inc1 = 1;
    inc2 = np0_;
    istart = k * np0_ * np1_;
    fftw(bwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
    // Transform second block along x: negative y indices
    inc1 = 1;
    inc2 = np0_;
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    fftw(bwplan0,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);

    // transform along y for all values of x
    ntrans = np0_;
    inc1 = np0_;
    inc2 = 1;
    istart = k * np0_ * np1_;
    fftw(bwplan1,ntrans,(FFTW_COMPLEX*)&val[istart],inc1,inc2,
                        (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP
  } // k
#elif defined(FFT_NOLIB) // USE_FFTW2
  // No library
  for ( int k = 0; k < np2_loc(); k++ )
  {
    // transform along x for non-zero vectors only
    // transform along x for y in [0,ntrans0_] and y in [np1-ntrans0_, np1-1]
    // transform along x for non-zero elements
    // Transform first block along x: positive y indices
    int ntrans = ntrans0_;
    int istart = k * np0_ * np1_;
    int length = np0_;
    int ainc   = 1;
    int ajmp   = np0_;
    double scale = 1.0;
    int idir = -1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // Transform second block along x: negative y indices
    istart = np0_ * ( (np1_-ntrans) + k * np1_ );
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );

    // transform along y for all values of x
    ntrans = np0_;
    istart = k * np0_ * np1_;
    length = np1_;
    ainc = np0_;
    ajmp = 1;
    cfftm (&val[istart],&val[istart],scale,ntrans,length,ainc,ajmp,idir );
  } // for k
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_fxy_inv.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fz(void)
{
  // transform zvec along z
#if TIMING
  tm_fz.start();
#endif

#if USE_ESSL_FFT
  int inc1 = 1, inc2 = np2_, ntrans = nvec_, isign = 1, initflag = 0;
  double scale = 1.0 / np012();

  if ( ntrans > 0 )
    dcft_(&initflag,&zvec_[0],&inc1,&inc2,&zvec_[0],&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);

#elif USE_FFTW2
#if _OPENMP
  const double fac = 1.0 / np012();
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    //#pragma omp task
    fftw_one(fwplan2,(FFTW_COMPLEX*)&zvec_[i*np2_],(FFTW_COMPLEX*)0);
    for ( int j = 0; j < np2_; j++ )
      zvec_[i*np2_+j] *= fac;
  }
  // int inc1=1;
  // zdscal(&len,&fac,&zvec_[0],&inc1);
#else

 /*
  * void fftw(fftw_plan plan, int howmany,
  *    FFTW_COMPLEX *in, int istride, int idist,
  *    FFTW_COMPLEX *out, int ostride, int odist);
  */
  int ntrans, inc1, inc2;

  ntrans = nvec_;
  inc1 = 1;
  inc2 = np2_;
  fftw(fwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
  int len = zvec_.size();
  double fac = 1.0 / np012();
  zdscal(&len,&fac,&zvec_[0],&inc1);
#endif
#elif USE_FFTW3

#if USE_FFTW3_THREADS
  fftw_execute_dft ( fwplan, (fftw_complex*)&zvec_[0],
                    (fftw_complex*)&zvec_[0]);
#else
  // do np2_ same for D_USE_1D or not
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_execute_dft ( fwplan, (fftw_complex*)&zvec_[i*np2_],
                      (fftw_complex*)&zvec_[i*np2_]);
  }
#endif
  // scale
  double fac = 1.0 / np012();
  int len = zvec_.size();
  int inc1 = 1;
  zdscal(&len,&fac,&zvec_[0],&inc1);
#elif defined(FFT_NOLIB)
  // No library
  /* Transform along z */
  int ntrans = nvec_;
  int length = np2_;
  int ainc   = 1;
  int ajmp   = np2_;
  double scale = 1.0 / np012();
  int idir = 1;
  cfftm ( &zvec_[0], &zvec_[0], scale, ntrans, length, ainc, ajmp, idir );
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_fz.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::fz_inv(void)
{
  // transform zvec along z
#if TIMING
  tm_fz_inv.start();
#endif

#if USE_ESSL_FFT
  int inc1 = 1, inc2 = np2_, ntrans = nvec_, isign = -1, initflag = 0;
  double scale = 1.0;

  if ( ntrans > 0 )
    dcft_(&initflag,&zvec_[0],&inc1,&inc2,&zvec_[0],&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
#elif USE_FFTW2
   /*
    * void fftw(fftw_plan plan, int howmany,
    *    FFTW_COMPLEX *in, int istride, int idist,
    *    FFTW_COMPLEX *out, int ostride, int odist);
    */
#if _OPENMP
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_one(bwplan2,(FFTW_COMPLEX*)&zvec_[i*np2_],(FFTW_COMPLEX*)0);
  }
#else
  int ntrans = nvec_, inc1 = 1, inc2 = np2_;
  fftw(bwplan2,ntrans,(FFTW_COMPLEX*)&zvec_[0],inc1,inc2,
                      (FFTW_COMPLEX*)0,0,0);
#endif // _OPENMP

#elif USE_FFTW3 // USE_FFTW2

#if USE_FFTW3_THREADS
  fftw_execute_dft ( bwplan, (fftw_complex*)&zvec_[0],
                     (fftw_complex*)&zvec_[0]);
#else
  #pragma omp parallel for
  for ( int i = 0; i < nvec_; i++ )
  {
    fftw_execute_dft ( bwplan, (fftw_complex*)&zvec_[i*np2_],
                       (fftw_complex*)&zvec_[i*np2_]);
  }
#endif // USE_FFTW3_THREADS

#elif defined(FFT_NOLIB) // USE_FFTW3
  // No library
  /* Transform along z */
  int ntrans = nvec_;
  int length = np2_;
  int ainc   = 1;
  int ajmp   = np2_;
  double scale = 1.0;
  int idir = -1;
  cfftm ( &zvec_[0], &zvec_[0], scale, ntrans, length, ainc, ajmp, idir );
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif // USE_FFTW3

#if TIMING
  tm_fz_inv.stop();
#endif
}

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::init_lib(void)
{
  // initialization of FFT libs
#if TIMING
  tm_init.start();
#endif

#if USE_ESSL_FFT
  complex<double> *p = 0;
#if USE_ESSL_2DFFT
  // use 2D FFT for x and y transforms and 1D FFT for z transforms
  naux1xy = 40000 + 2.28 * (np0_+np1_);
  aux1xyf.resize(naux1xy);
  aux1xyb.resize(naux1xy);
  int r = max(np0_,np1_);
  int s = min(64,min(np0_,np1_));
  naux2 = 20000 + (2*r+256)*(s+2.28);

  naux1z = 20000 + 2.28 * np2_;
  aux1zf.resize(naux1z);
  aux1zb.resize(naux1z);

  int ntrans2 = nvec_;
  int naux2z = 20000 + 2.28 * np2_ + (256 + 2*np2_)*min(64,ntrans2);
  naux2 = max( naux2, naux2z );
  aux2.resize(naux2);

  double scale = 1.0;

  // initialize xy transforms
  int initflag = 1, inc1, inc2, isign = -1;
  inc1 = 1; inc2 = np0_;
  dcft2_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&np1_,
         &isign,&scale,&aux1xyb[0],&naux1xy,&aux2[0],&naux2);
  isign = 1;
  dcft2_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&np1_,
         &isign,&scale,&aux1xyf[0],&naux1xy,&aux2[0],&naux2);

  // initialize z transforms
  int ntrans = nvec_;
  if ( ntrans > 0 )
  {
    inc1 = 1; inc2 = np2_;
    isign = -1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
    isign = 1; scale = 1.0 / np012();
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);
  }
#else // USE_ESSL_2DFFT

  naux1x = (int) (20000 + 2.28 * np0_);
  naux1y = (int) (20000 + 2.28 * np1_);
  naux1z = (int) (20000 + 2.28 * np2_);
  aux1xf.resize(naux1x);
  aux1yf.resize(naux1y);
  aux1zf.resize(naux1z);
  aux1xb.resize(naux1x);
  aux1yb.resize(naux1y);
  aux1zb.resize(naux1z);

  int naux2x = (int) (20000 + 2.28 * np0_ + (256 + 2*np0_)*min(64,ntrans0_));
  naux2 = naux2x;
  int naux2y = (int) (20000 + 2.28 * np1_ + (256 + 2*np1_)*min(64,ntrans1_));
  naux2 = max( naux2, naux2y );
  int naux2z = (int) (20000 + 2.28 * np2_ + (256 + 2*np2_)*min(64,ntrans2_));
  naux2 = max( naux2, naux2z );
  aux2.resize(naux2);

  // initialize x, y and z transforms

  int initflag = 1, inc1, inc2, ntrans, isign;
  double scale = 1.0;

  // x transforms
  inc1 = 1; inc2 = np0_; ntrans = ntrans0_;
  if ( ntrans > 0 )
  {
    isign = -1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&ntrans,
          &isign,&scale,&aux1xb[0],&naux1x,&aux2[0],&naux2);
    isign = 1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np0_,&ntrans,
          &isign,&scale,&aux1xf[0],&naux1x,&aux2[0],&naux2);
  }

  // y transforms
  inc1 = np0_; inc2 = 1; ntrans = ntrans1_;
  if ( ntrans > 0 )
  {
    isign = -1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np1_,&ntrans,
          &isign,&scale,&aux1yb[0],&naux1y,&aux2[0],&naux2);
    isign = 1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np1_,&ntrans,
          &isign,&scale,&aux1yf[0],&naux1y,&aux2[0],&naux2);
  }

  // z transforms
  inc1 = 1; inc2 = np2_; ntrans = ntrans2_;
  if ( ntrans > 0 )
  {
    isign = -1;
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zb[0],&naux1z,&aux2[0],&naux2);
    isign = 1; scale = 1.0 / np012();
    dcft_(&initflag,p,&inc1,&inc2,p,&inc1,&inc2,&np2_,&ntrans,
          &isign,&scale,&aux1zf[0],&naux1z,&aux2[0],&naux2);
  }

#endif // USE_ESSL_2DFFT

#elif USE_FFTW2

  fwplan0 = fftw_create_plan(np0_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  fwplan1 = fftw_create_plan(np1_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  fwplan2 = fftw_create_plan(np2_,FFTW_FORWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan0 = fftw_create_plan(np0_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan1 = fftw_create_plan(np1_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);
  bwplan2 = fftw_create_plan(np2_,FFTW_BACKWARD,FFTW_ALGO|FFTW_IN_PLACE);

#elif USE_FFTW3
  vector<complex<double> > aux(np0_*np1_);
#if defined(USE_FFTW3MKL) && !defined(USE_FFTW3_THREADS) && _OPENMP
  fftw3_mkl.number_of_user_threads = omp_get_max_threads();
#endif

#if USE_FFTW3_THREADS
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  vector<complex<double> > aux1(np0_*np1_*np2_loc());

  // xy
  int rank = 2;
  int n[] = {np1_,np0_};
  int howmany = np2_loc();
  //int howmany = 1;
  int idist = np0_*np1_, odist = np0_*np1_;
  int istride = 1, ostride = 1; /* array is contiguous in memory */
  int *inembed = n, *onembed = n;

  fwplan2d = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux1[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux1[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplan2d = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux1[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux1[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);

  // z
  rank = 1;
  int nz[] = {np2_};
  howmany = nvec_;
  idist = np2_, odist = np2_;
  istride = 1, ostride = 1; /* array is contiguous in memory */
  inembed = nz, onembed = nz;

  fwplan = fftw_plan_many_dft(rank, nz, howmany, (fftw_complex*)&zvec_[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&zvec_[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplan = fftw_plan_many_dft(rank, nz, howmany, (fftw_complex*)&zvec_[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&zvec_[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);


#else // USE_FFTW3_THREADS
#if USE_FFTW3_2D
  // row major in FFTW3 2d plans
  fwplan2d = fftw_plan_dft_2d ( np1_, np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplan2d = fftw_plan_dft_2d ( np1_, np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );
#else // USE_FFTW3_2D
  // FFTW3 1D
  fwplanx = fftw_plan_dft_1d ( np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplanx = fftw_plan_dft_1d ( np0_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );

#if FFTW_TRANSPOSE
  fwplany = fftw_plan_dft_1d ( np1_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), -1,
                                FFTW_ALGO );
  bwplany = fftw_plan_dft_1d ( np1_, (fftw_complex*)(&aux[0]),
                                (fftw_complex*)(&aux[0]), 1,
                                FFTW_ALGO );

#else // FFTW_TRANSPOSE
  // strided FFT
  int rank = 1;
  int n[] = {np1_};
  int howmany = 1;
  int idist = 1, odist = 1;
  int istride = np0_, ostride = np0_; /* array is contiguous in memory */
  int *inembed = n, *onembed = n;

  fwplany = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux[0], onembed,
                                    ostride, odist, -1, FFTW_ALGO);
  bwplany = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)&aux[0],
                                    inembed, istride, idist,
                                    (fftw_complex*)&aux[0], onembed,
                                    ostride, odist, 1, FFTW_ALGO);
#endif // FFTW_TRANSPOSE
#endif // USE_FFTW3_2D
  // do z using 1d plans
  fwplan = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zvec_[0]),
                                (fftw_complex*)(&zvec_[0]), -1,
                                FFTW_ALGO );
  bwplan = fftw_plan_dft_1d ( np2_, (fftw_complex*)(&zvec_[0]),
                                (fftw_complex*)(&zvec_[0]), 1,
                                FFTW_ALGO );
#endif //USE_FFTW3_THREADS

#elif FFT_NOLIB // USE_FFTW3
  /* no library */
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if TIMING
  tm_init.stop();
#endif
}

#if defined(FFT_NOLIB)

////////////////////////////////////////////////////////////////////////////////
//
//     /* no library: use cfftm function */
//
//
//     /* Transform along x */
//     int i2;
//     int ntrans = np1*np2;
//     int length = np0;
//     int ainc   = 1;
//     int ajmp   = np0i;
//     int idir;
//     double scale;
//     if ( dir == R_TO_K )
//     {
//       idir = -1;
//       scale = 1.0 / ((double) np0*np1*np2);
//     }
//     else
//     {
//       idir = 1;
//       scale = 1.0;
//     }
//
//     cfftm ( &val[0], &val[0], scale, ntrans, length, ainc, ajmp, idir );
//
//     /* Transform along y */
//     for ( i2 = 0; i2 < np2; i2++ )
//     {
//       int ist = i2 * np0i * np1;
//       ntrans = np0;
//       length = np1;
//       ainc   = np0i;
//       ajmp   = 1;
//       scale = 1.0;
//       cfftm ( &val[ist], &val[ist], scale, ntrans, length, ainc, ajmp, idir );
//     }
//
//     /* Transform along z */
//     ntrans = np0i*np1;
//     length = np2;
//     ainc   = np0i*np1;
//     ajmp   = 1;
//     scale = 1.0;
//     cfftm ( &val[0], &val[0], scale, ntrans, length, ainc, ajmp, idir );
//
////////////////////////////////////////////////////////////////////////////////

/* define multiple FFT function here */

void cfftm ( complex<double> *ain, complex<double> *aout, double scale,
  int ntrans, int length,
  int ainc, int ajmp, int idir )
/*
 *  cfftm: multiple one-dimensional complex FFT
 *
 *  ain     complex array (input)
 *  aout    complex array (output)
 *  scale   global scaling factor
 *  ntrans  number of transforms
 *  length  length of each transform (in complex numbers)
 *  ainc    distance between elements within a transform (in complex numbers)
 *  ajmp    distance between first elements of transforms (in complex numbers)
 *  idir    direction of transform
 */

{
  void cfft ( int idir, complex<double> *z1, complex<double> *z2, int n,
    int *inzee );
  vector<complex<double> > tmpa(length), tmpb(length);
  for ( int it = 0; it < ntrans; it++ )
  {
    int ibase = it * ajmp;
    for ( int i = 0; i < length; i++ )
    {
      tmpa[i] = ain[ibase+i*ainc];
    }
    int inzee = 1;
    cfft ( idir, &tmpa[0], &tmpb[0], length, &inzee );
    if ( inzee == 1 )
      for ( int i = 0; i < length; i++ )
      {
        aout[ibase+i*ainc] = tmpa[i];
      }
    else
      for ( int i = 0; i < length; i++ )
      {
        aout[ibase+i*ainc] = tmpb[i];
      }
    for ( int i = 0; i < length; i++ )
    {
      aout[ibase+i*ainc] *= scale;
    }
  }
}

/*******************************************************************************
 *
 *  Complex FFT
 *  C version
 *
 *  From: C++ Language System Release 3.0 Library Manual
 *  Transcription from 'FFT as Nested Multiplication, with a twist'
 *  C. de Boor, SIAM Sci. Stat. Comput., Vol 1, No 1, March 1980
 *
 *  Adapted to C  17 Feb 1993, 9 Dec 1993
 *
 ******************************************************************************/

#include <math.h>

#define NEXTMX 12

void cfft ( int idir, complex<double> *z1, complex<double> *z2, int n,
  int *inzee )
{
  // Compute the discrete Fourier transform of z1 (or z2) in
  // the Cooley-Tukey way, but with a twist.
  // z1[before], z2[before]
  // *inzee == 1 means input in z1; *inzee == 2 means input in z2

  void fftstp ( int idir, complex<double> *zin, int after,
                int now, int before, complex<double> *zout );
  static int prime[NEXTMX] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
  int before = n;
  int after = 1;
  int next = 0;
  int now;


  do
  {
    int np = prime[next];
    if ( (before/np)*np < before )
    {
      if ( ++next < NEXTMX ) continue;
      now = before;
      before = 1;
    }
    else
    {
      now = np;
      before /= np;
    }
    if ( *inzee == 1 )
      fftstp ( idir, z1, after, now, before, z2 );
    else
      fftstp ( idir, z2, after, now, before, z1 );
    *inzee = 3 - *inzee;
    after *= now;
  } while ( before > 1 );

}

void fftstp ( int idir, complex<double> *zin, int after,
              int now, int before, complex<double> *zout )
{

  static const double twopi = 2 * 3.141592653589793;
  double angle;
  complex<double> omega;
  complex<double> arg,value;
  int ia,ib,j,k,l,in;

  angle = twopi/(now*after);
  omega =  complex<double>(cos ( angle ),-idir * sin ( angle ));
  arg = 1.0;
  for ( j = 0; j < now; j++ )
  {
    for ( ia = 0; ia < after; ia++ )
    {
      for ( ib = 0; ib < before; ib++ )
      {
        /* value = zin(ia,ib,now) */
        k = ia + ib*after + (now-1)*before*after;
        value = zin[k];
        for ( in = now-2; in >= 0; in-- )
        {
          /* value = value*arg + zin(ia,ib,in) */
          /* zin(ia,ib,in) = zin[ia + ib*after + in*before*after]; */
          l = ia + ib*after + in*before*after;
          value = value * arg + zin[l];
        }
        /* zout(ia,j,ib) = value */
        /* zout[ia + j*after + ib*now*after] = value; */
        l = ia + j*after + ib*now*after;
        zout[l] = value;
      }

      /* arg *= omega; */
      arg *= omega;

    }
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
void FourierTransform::reset_timers(void)
{
#if TIMING
  tm_fwd.reset();
  tm_bwd.reset();
  tm_map_fwd.reset();
  tm_map_fwd.reset();
  tm_trans_fwd.reset();
  tm_trans_bwd.reset();
  tm_fxy.reset();
  tm_fxy_inv.reset();
  tm_fz.reset();
  tm_fz_inv.reset();
  tm_init.reset();
#endif
}
