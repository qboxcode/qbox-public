////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// D3tensor.h
//
// double 3x3 tensor
//
////////////////////////////////////////////////////////////////////////////////

#ifndef D3TENSOR_H
#define D3TENSOR_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include "blas.h"
#include "D3vector.h"

using namespace std;

class D3tensor
{
  public:

  double a_[9];

  double* a(void) { return &a_[0]; }
  const double* a(void) const { return &a_[0]; }

  explicit D3tensor(void) { clear(); }

  explicit D3tensor(double xx, double yy, double zz)
  { a_[0]=xx; a_[4]=yy; a_[8]=zz; }

  explicit D3tensor(double xx, double yy, double zz,
                    double xy, double yz, double xz, char& uplo)
  {
    a_[0] = xx;
    a_[4] = yy;
    a_[8] = zz;

    if ( uplo == 'l' )
    {
      a_[1] = xy;
      a_[2] = xz;
      a_[5] = yz;
    }
    else if ( uplo == 'u' )
    {
      a_[3] = xy;
      a_[6] = xz;
      a_[7] = yz;
    }
    else if ( uplo == 's' )
    {
      a_[1] = xy;
      a_[2] = xz;
      a_[5] = yz;
      a_[3] = xy;
      a_[6] = xz;
      a_[7] = yz;
    }
    else
      assert(false);
  }

  explicit D3tensor(const D3vector& diag, const D3vector& offdiag)
  {
    a_[0] = diag[0];
    a_[4] = diag[1];
    a_[8] = diag[2];
    a_[1] = offdiag[0];
    a_[5] = offdiag[1];
    a_[2] = offdiag[2];
    a_[3] = offdiag[0];
    a_[7] = offdiag[1];
    a_[6] = offdiag[2];
  }

  explicit D3tensor(const double* a)
  {
    for ( int i = 0; i < 9; i++ ) a_[i] = a[i];
  }

  double& operator[](int i)
  {
    assert(i>=0 && i < 9);
    return a_[i];
  }

  double operator[](int i) const
  {
    assert(i>=0 && i < 9);
    return a_[i];
  }

  void setdiag(int i, double b)
  {
    assert(i>=0 && i<3);
    a_[i*4] = b;
  }

  void setdiag(const D3vector& b)
  {
    for ( int i = 0; i < 3; i++ )
      a_[i*4] = b[i];
  }

  void setoffdiag(int i, double b)
  {
    assert(i>=0 && i<3);
    if ( i == 0 )
    {
      a_[1] = b;
      a_[3] = b;
    }
    else if ( i == 2 )
    {
      a_[2] = b;
      a_[6] = b;
    }
    else
    {
      a_[5] = b;
      a_[7] = b;
    }
  }

  void setoffdiag(const D3vector& b)
  {
    a_[1] = b[0];
    a_[3] = b[0];
    a_[5] = b[1];
    a_[7] = b[1];
    a_[2] = b[2];
    a_[6] = b[2];
  }

  bool operator==(const D3tensor &rhs) const
  {
    bool eq = true;
    for ( int i = 0; i < 9; i++ )
    {
      if ( rhs[i] != a_[i] )
      {
        eq &= false;
      }
    }
    return eq;
  }

  bool operator!=(const D3tensor &rhs) const
  {
    bool neq = false;
    for ( int i = 0; i < 9; i++ )
    {
      if ( rhs[i] != a_[i] )
      {
        neq |= true;
      }
    }
    return neq;
  }

  D3tensor& operator+=(const D3tensor& rhs)
  {
    for ( int i = 0; i < 9; i++ )
      a_[i] += rhs[i];
    return *this;
  }

  D3tensor& operator-=(const D3tensor& rhs)
  {
    for ( int i = 0; i < 9; i++ )
      a_[i] -= rhs[i];
    return *this;
  }

  D3tensor& operator*=(const double& rhs)
  {
    for ( int i = 0; i < 9; i++ )
      a_[i] *= rhs;
    return *this;
  }

  D3tensor& operator/=(const double& rhs)
  {
    for ( int i = 0; i < 9; i++ )
      a_[i] /= rhs;
    return *this;
  }

  friend const D3tensor operator+(const D3tensor& lhs, const D3tensor& rhs)
  {
    return D3tensor(lhs) += rhs;
  }

  friend const D3tensor operator-(const D3tensor& a, const D3tensor& b)
  {
    return D3tensor(a) -= b;
  }

  friend D3tensor operator*(double a, const D3tensor& b)
  {
    return D3tensor(b) *= a;
  }

  friend D3tensor operator*(const D3tensor& a, double b)
  {
    return D3tensor(a) *= b;
  }

  friend D3tensor operator/(const D3tensor& a, double b)
  {
    return D3tensor(a) /= b;
  }

  friend D3tensor operator*(D3tensor& a, D3tensor& b)
  {
    D3tensor c;
    int ithree = 3;
    double one = 1.0, zero = 0.0;
    char t = 'n';
    dgemm ( &t, &t, &ithree, &ithree, &ithree, &one, &a[0], &ithree,
            &b[0], &ithree, &zero, &c[0], &ithree );
    return c;
  }

  friend D3tensor operator-(const D3tensor& a) // unary minus
  {
    return D3tensor()-a;
  }

  double norm2(const D3tensor& a) const
  {
    return a[0]*a[0] + a[1]*a[1] + a[2]*a[2] +
           a[3]*a[3] + a[4]*a[4] + a[5]*a[5] +
           a[6]*a[6] + a[7]*a[7] + a[8]*a[8];
  }

  double norm(const D3tensor& a) const
  {
    return sqrt(norm2(a));
  }

  double trace(void) const
  {
    return a_[0]+a_[4]+a_[8];
  }

  void traceless(void)
  {
    double b = trace() / 3;
    a_[0] -= b;
    a_[4] -= b;
    a_[8] -= b;
  }

  void clear(void)
  {
    for ( int i = 0; i < 9; i++ )
      a_[i] = 0.0;
  }

  void identity(void)
  {
    clear();
    a_[0] = 1.0;
    a_[4] = 1.0;
    a_[8] = 1.0;
  }

  friend std::ostream& operator<<(std::ostream& os, const D3tensor& rhs)
  {
    const double * const v  = rhs.a();
    os.setf(ios::fixed,ios::floatfield);
    os.setf(ios::right,ios::adjustfield);
    os.precision(8);
    os << setw(14) << v[0] << " " << setw(14) << v[3] << " " << setw(14) << v[6]
       << "\n"
       << setw(14) << v[1] << " " << setw(14) << v[4] << " " << setw(14) << v[7]
       << "\n"
       << setw(14) << v[2] << " " << setw(14) << v[5] << " " << setw(14) << v[8]
       << "\n";
    return os;
  }

  void syev(char uplo, D3vector& eigval, D3tensor& eigvec)
  {
    double w[3];

    int info;
    char jobz = 'V';
    int lwork=-1;
    double tmplwork;
    int n = 3;

    dsyev(&jobz, &uplo, &n, eigvec.a(), &n, &w[0], &tmplwork, &lwork, &info);

    lwork = (int) tmplwork + 1;
    double* work = new double[lwork];

    eigvec = *this;
    dsyev(&jobz, &uplo, &n, eigvec.a(), &n, &w[0], work, &lwork, &info);
    delete[] work;

    eigval = D3vector(&w[0]);
  }
};
#endif
