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
// FourierTransform.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FOURIERTRANSFORM_H
#define FOURIERTRANSFORM_H

#include <complex>
#include <vector>

#if !( defined(USE_FFTW2) || defined(USE_FFTW3) || defined(USE_ESSL_FFT) || defined(FFT_NOLIB) )
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

#if defined(USE_FFTW2) && defined(USE_FFTW3)
#error "Cannot define USE_FFTW2 and USE_FFTW3"
#endif

#if USE_FFTW2
#if USE_DFFTW
#include "dfftw.h"
#else
#include "fftw.h"
#endif
#endif

#if USE_FFTW3
#include "fftw3.h"
#if USE_FFTW3MKL
#include "fftw3_mkl.h"
#endif
#endif

#include "Timer.h"
#include "BasisMapping.h"

class Basis;

class FourierTransform
{
  private:

  const Basis& basis_;
  const BasisMapping bm_;

  const int np0_,np1_,np2_;
  const int nvec_;

  int ntrans0_,ntrans1_,ntrans2_;

  std::vector<std::complex<double> > zvec_;
  void init_lib(void);

#if USE_ESSL_FFT
#if USE_ESSL_2DFFT
  std::vector<double> aux1xyf,aux1zf;
  std::vector<double> aux1xyb,aux1zb;
  std::vector<double> aux2;
  int naux1xy,naux1z,naux2;
#else
  std::vector<double> aux1xf, aux1yf, aux1zf;
  std::vector<double> aux1xb, aux1yb, aux1zb;
  std::vector<double> aux2;
  int naux1x,naux1y,naux1z,naux2;
#endif
#elif USE_FFTW2
  fftw_plan fwplan0,fwplan1,fwplan2,bwplan0,bwplan1,bwplan2;
#elif USE_FFTW3
  //plans for np2_
  fftw_plan fwplan, bwplan;
#if defined(USE_FFTW3_2D) || defined(USE_FFTW3_THREADS)
  fftw_plan fwplan2d, bwplan2d;
#else
  fftw_plan fwplanx, fwplany, bwplanx, bwplany;
#endif
#elif defined(FFT_NOLIB)
  // no library
#else
#error "Must define USE_FFTW2, USE_FFTW3, USE_ESSL_FFT or FFT_NOLIB"
#endif

  void fwd(std::complex<double>* val);
  void bwd(std::complex<double>* val);

  public:

  FourierTransform (const Basis &basis, int np0, int np1, int np2);
  ~FourierTransform ();

  // backward: Fourier synthesis, compute real-space function
  // forward:  Fourier analysis, compute Fourier coefficients
  // forward transform includes scaling by 1/np012
  // single transforms: c -> f, f -> c
  void backward (const std::complex<double>* c, std::complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f, std::complex<double>* c);

  // double transforms: c1 + i*c2 -> f, f -> c1 + i*c2
  void backward (const std::complex<double>* c1,
                 const std::complex<double>* c2, std::complex<double>* f);
  // Note: forward transforms overwrite the array f
  void forward(std::complex<double>* f,
               std::complex<double>* c1, std::complex<double>* c2);

  int np0() const { return np0_; }
  int np1() const { return np1_; }
  int np2() const { return np2_; }
  int np2_loc(void) const { return bm_.np2_loc(); }
  int np2_loc(int iproc) const { return bm_.np2_loc(iproc); }
  int np2_first(void) const { return bm_.np2_first(); }
  int np2_first(int iproc) const { return bm_.np2_first(iproc); }
  long int np012() const { return ((long int)np0_) * np1_ * np2_; }
  int np012loc(int iproc) const { return np0_ * np1_ * np2_loc(iproc); }
  int np012loc(void) const { return np0_ * np1_ * np2_loc(); }
  int index(int i, int j, int k) const
  { return i + np0_ * ( j +  np1_ * k ); }
  int i(int ind) const { return ind % np0_; }
  int j(int ind) const { return (ind / np0_) % np1_; }
  int k(int ind) const { return (ind / np0_) / np1_ + np2_first(); }

  void reset_timers(void);
  Timer tm_f_map, tm_f_fft, tm_f_pack, tm_f_mpi, tm_f_zero, tm_f_unpack,
        tm_b_map, tm_b_fft, tm_b_pack, tm_b_mpi, tm_b_zero, tm_b_unpack,
        tm_f_xy, tm_f_z, tm_f_x, tm_f_y,
        tm_b_xy, tm_b_z, tm_b_x, tm_b_y,
        tm_init, tm_b_com, tm_f_com;

};
#endif
