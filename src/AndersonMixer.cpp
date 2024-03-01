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
// AndersonMixer.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "AndersonMixer.h"
#include "blas.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AndersonMixer::AndersonMixer(const int m, const int nmax,
  bool distributed) : m_(m), nmax_(nmax), diag_(true),
  eig_ratio_(0.0), distributed_(distributed)
{
  assert( nmax >= 0 );
  x_.resize(nmax_+1);
  f_.resize(nmax_+1);
  for ( int n = 0; n < nmax_+1; n++ )
  {
    x_[n].resize(m_);
    f_[n].resize(m_);
  }
  restart();
}

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::restart(void)
{
  n_ = -1;
  k_ = -1;
}

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::update(double* x, double* f, double* xbar, double* fbar)
{
  // update:
  // input: x, f
  // output: xbar, fbar
  //
  // Computes the pair (xbar,fbar) using pairs (x,f) used
  // in previous updates, according to the Anderson algorithm.

  // increment index of current vector
  k_ = ( k_ + 1 ) % ( nmax_ + 1 );
  // increment current number of vectors
  if ( n_ < nmax_ ) n_++;

  // save vectors

  for ( int i = 0; i < m_; i++ )
  {
    x_[k_][i] = x[i];
    f_[k_][i] = f[i];
  }

  valarray<double> a,atmp;
  valarray<double> b,btmp;
  valarray<double> theta;
  if ( n_ > 0 )
  {
    // compute matrix A = F^T F and rhs b = F^T f
    // compute the lower part of A only (i>=j)
    a.resize(n_*n_);
    atmp.resize(n_*n_);
    b.resize(n_);
    btmp.resize(n_);
    theta.resize(n_);
    for ( int i = 0; i < n_; i++ )
    {
      const int kmi = ( k_ - i + nmax_ ) % ( nmax_ + 1 );
      assert(kmi>=0);
      assert(kmi<nmax_+1);
      for ( int j = 0; j <= i; j++ )
      {
        const int kmj = ( k_ - j + nmax_ ) % ( nmax_ + 1 );
        assert(kmj>=0);
        assert(kmj<nmax_+1);
        double sum = 0.0;
        for ( int l = 0; l < m_; l++ )
          sum += (f_[k_][l] - f_[kmi][l]) * (f_[k_][l] - f_[kmj][l]);
        a[i+j*n_] = sum;
      }
      double bsum = 0.0;
      for ( int l = 0; l < m_; l++ )
        bsum += ( f_[k_][l] - f_[kmi][l] ) * f_[k_][l];
      b[i] = bsum;
    }

    if ( distributed_ )
    {
      MPI_Allreduce(&a[0],&atmp[0],n_*n_,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
      a = atmp;
      MPI_Allreduce(&b[0],&btmp[0],n_,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
      b = btmp;
    }

    // solve the linear system a * theta = b
    // solve on task 0 and bcast result
    if ( !distributed_ || MPIdata::rank() == 0 )
    {
      if ( diag_ )
      {
        // solve the linear system using eigenvalues and eigenvectors
        // compute eigenvalues of a
        char jobz = 'V';
        char uplo = 'L';
        valarray<double> w(n_);
        int lwork = 3*n_;
        valarray<double> work(lwork);
        int info;
        dsyev(&jobz,&uplo,&n_,&a[0],&n_,&w[0],&work[0],&lwork,&info);

#ifdef DEBUG
        cout << "AndersonMixer: eigenvalues: ";
        for ( int i = 0; i < n_; i++ )
          cout << w[i] << "  ";
        cout << endl;
#endif
        if ( info != 0 )
        {
          cerr << " AndersonMixer: Error in dsyev" << endl;
          cerr << " info = " << info << endl;
          exit(1);
        }

        // solve for theta
        // inverse(a) = sum_k (1/w[k]) * u_k * u_k^T
        // the eigenvector u_k is the k_th column of the matrix a
        // theta = inverse(a) * b = sum_k (1/w[k]) * u_k * ( u_k^T * b )
        // set all elements of theta to 0.0
        theta = 0.0;
        for ( int k = 0; k < n_; k++ )
        {
          // use u_k only if eigenvalue w[k] is large enough compared to the
          // largest eigenvalue
          if ( w[k] > eig_ratio_ * w[n_-1] )
          {
            // compute dot product u_k^T * b
            double uktb = 0.0;
            for ( int j = 0; j < n_; j++ )
              uktb += a[j+k*n_] * b[j];
            for ( int i = 0; i < n_; i++ )
              theta[i] += ( uktb / w[k] ) * a[i+k*n_];
          }
        }
      }
      else
      {
        // solve the linear system directly

        // Tikhonov regularization parameter
        // adjust the parameter until the norm of theta is < 1.0
        double tikhonov_parameter = 1.e-12;
        bool norm_ok = false;
        valarray<double> asave(a);
        valarray<double> bsave(b);
        int iter = 0;
        const int maxiter = 100;
        while ( !norm_ok && iter < maxiter )
        {
          a = asave;
          b = bsave;
          for ( int i = 0; i < n_; i++ )
            a[i+i*n_] += tikhonov_parameter;

          char uplo = 'L';
          int nrhs = 1;
          valarray<int> ipiv(n_);
          valarray<double> work(n_);
          int info;
          dsysv(&uplo,&n_,&nrhs,&a[0],&n_,&ipiv[0],
                &b[0],&n_,&work[0],&n_,&info);
          if ( info != 0 )
          {
            cerr << " AndersonMixer: Error in dsysv" << endl;
            cerr << " info = " << info << endl;
            exit(1);
          }
          // the vector b now contains the solution
          theta = b;

          // check condition on the norm of theta
          norm_ok = true;
#ifdef ANDERSON_SIMPLEX
          // unit simplex criterion
          double theta_sum = 0.0;
          for ( int i = 0; i < theta.size(); i++ )
          {
            theta_sum += theta[i];
            norm_ok &= theta[i] >= 0.0;
          }
          norm_ok &= fabs(theta_sum) <= 1.0;
#endif
#ifdef ANDERSON_INF_NORM
          // infinity norm criterion
          for ( int i = 0; i < theta.size(); i++ )
            norm_ok &= fabs(theta[i]) <  3.0;
#endif
#ifdef ANDERSON_2_NORM
          // 2-norm criterion
          double theta_sum = 0.0;
          for ( int i = 0; i < theta.size(); i++ )
          {
            theta_sum += theta[i] * theta[i];
          }
          norm_ok = theta_sum <= 1.0;
#endif
#ifdef DEBUG
          cout << " tp = " << tikhonov_parameter
               << " AndersonMixer: theta = ";
          for ( int i = 0; i < theta.size(); i++ )
            cout << theta[i] << " ";
          cout << endl;
#endif

          tikhonov_parameter *= 2.0;
          iter++;
        }
      }

#ifdef DEBUG
      cout << " AndersonMixer: theta = ";
      for ( int i = 0; i < theta.size(); i++ )
        cout << theta[i] << " ";
      cout << endl;
#endif

    }

    // broadcast theta from task 0
    if ( distributed_ )
    {
      MPI_Bcast(&theta[0],n_,MPI_DOUBLE,0,MPIdata::comm());
    }

  } // if n_ > 0

  // fbar = f[k] + sum_i theta_i * ( f[k] - f[kmi] )
  // xbar = x[k] + sum_i theta_i * ( x[k] - x[kmi] )
  for ( int l = 0; l < m_; l++ )
  {
    fbar[l] = f_[k_][l];
    xbar[l] = x_[k_][l];
  }
  for ( int i = 0; i < n_; i++ )
  {
    const int kmi = ( k_ - i + nmax_ ) % ( nmax_ + 1 );
    assert(kmi>=0);
    assert(kmi<nmax_+1);
    for ( int l = 0; l < m_; l++ )
    {
      fbar[l] -= theta[i] * ( f_[k_][l] - f_[kmi][l] );
      xbar[l] -= theta[i] * ( x_[k_][l] - x_[kmi][l] );
    }
  }
}
