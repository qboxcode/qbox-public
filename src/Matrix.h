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
// Matrix.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#define MATRIX_DEF_BLOCK_SIZE 64

#include "Context.h"

#include <valarray>
#include <complex>
#include <cstring> // memcpy

int numroc0(int n, int nb, int iproc, int nprocs);

class ComplexMatrix;

class DoubleMatrix
{
  private:

    const Context& ctxt_;
    int ictxt_;
    int lld_;  // leading dimension of local matrix
    int m_, n_;   // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_;   // size of local array
    int size_; // size of local array;
    int nprow_, npcol_; // number of process rows and cols in context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block
    bool reference_; // object was created using the copy constructor
    double* val;

  public:

    double* valptr(int i=0) { return &val[i]; }
    const double* cvalptr(int i=0) const { return &val[i]; }
    double& operator[] (int i) { return val[i]; }
    const double& operator[] (int i) const { return val[i]; }
    const Context& context(void) const { return ctxt_; }
    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int mloc(int irow) const; // size of local array on row irow
    int nloc(void) const { return nloc_; } // size of local array
    int nloc(int icol) const; // size of local array on column icol
    int size(void) const { return size_; } // size of local array
    int localsize(void) const { return mloc_*nloc_; } // local size of val
    double memsize(void) const { return (double)m_*(double)n_*sizeof(double); }
    double localmemsize(void) const
    { return (double) mloc_ * (double) nloc_ * sizeof(double); }
    const int* desc(void) const { return &desc_[0]; }

    // local block size of block (l,m)
    int mbs(int l) const
    {
      // if l is the last block and this process holds an incomplete
      // block, then the size is m_%mb_. Otherwise, it is mb_.
      return ( (l==mblocks_-1) && ( m_incomplete_ ) ) ? m_%mb_ : mb_;
    }
    int nbs(int m) const
    {
      // if m is the last block and this process holds an incomplete
      // block, then the size is n_%nb_. Otherwise, it is nb_.
      return ( (m==nblocks_-1) && ( n_incomplete_ ) ) ? n_%nb_ : nb_;
    }

    // number of local blocks
    int mblocks(void) const { return mblocks_; }
    int nblocks(void) const { return nblocks_; }

    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)
    int l(int i) const { return i/(nprow_*mb_); }
    int x(int i) const { return i % mb_; }
    int pr(int i) const { return (i/mb_) % nprow_; }

    int m(int j) const { return j/(npcol_*nb_); }
    int y(int j) const { return j % nb_; }
    int pc(int j) const { return (j/nb_) % npcol_; }

    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(int l, int x) const { return (l * nprow_ + myrow_) * mb_ + x; }
    int j(int m, int y) const { return (m * npcol_ + mycol_) * nb_ + y; }

    int iglobal(int ilocal) const
    { return mb_*(nprow_*(ilocal/mb_)+myrow_)+ilocal%mb_; }
    int jglobal(int jlocal) const
    { return nb_*(npcol_*(jlocal/nb_)+mycol_)+jlocal%nb_; }

    int iglobal(int iprow, int ilocal) const
    { return mb_*(nprow_*(ilocal/mb_)+iprow)+ilocal%mb_; }
    int jglobal(int jpcol, int jlocal) const
    { return nb_*(npcol_*(jlocal/nb_)+jpcol)+jlocal%nb_; }


    // store element a(ii,jj) (where ii,jj are global indices)
    // in array val:
    //        int iii = l(ii) * mb_ + x(ii);
    //        int jjj = m(jj) * nb_ + y(jj);
    //        val[iii+mloc_*jjj] = a_ij;

    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }

    void init_size(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE,
      int nb = MATRIX_DEF_BLOCK_SIZE);

    void resize(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE,
      int nb = MATRIX_DEF_BLOCK_SIZE)
    {
      const int old_size = size_;
      init_size(m,n,mb,nb);
      if ( size_ == old_size ) return;
      delete[] val;
      val = new double[size_];
      clear();
    }

    void print(std::ostream& os) const;

    explicit DoubleMatrix(const Context& ctxt) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), mloc_(0), nloc_(0),
        size_(0), reference_(false), val(0) {}

    // Construct a DoubleMatrix of dimensions m,n
    explicit DoubleMatrix(const Context& ctxt, int m, int n,
        int mb=MATRIX_DEF_BLOCK_SIZE,
        int nb=MATRIX_DEF_BLOCK_SIZE) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0)
    {
      resize(m,n,mb,nb);
    }

    // copy constructor: create a separate copy of rhs
    explicit DoubleMatrix(const DoubleMatrix& rhs) : ctxt_(rhs.context()),
      reference_(false)
    {
      init_size(rhs.m(),rhs.n(),rhs.mb(),rhs.nb());
      val = new double[size_];
      memcpy(val, rhs.val, size_*sizeof(double));
    }

    // reference constructor create a proxy for a ComplexMatrix rhs
    explicit DoubleMatrix(ComplexMatrix& rhs);
    explicit DoubleMatrix(const ComplexMatrix& rhs);

    ~DoubleMatrix(void)
    {
      if ( !reference_ ) delete[] val;
    }

    DoubleMatrix& operator=(const DoubleMatrix& a);

    DoubleMatrix& operator+=(const DoubleMatrix& a);
    DoubleMatrix& operator-=(const DoubleMatrix& a);

    DoubleMatrix& operator*=(double a);

    void matgather(double *a, int lda) const;

    void initdiag(const double* const a); // initialize diagonal from a[]
    void init(const double* const a, int lda);
    void set(char uplo, double x);
    void identity(void);
    void clear(void);

    double dot(const DoubleMatrix &a) const;
    void axpy(double alpha, const DoubleMatrix &a);
    void scal(double alpha);

    double nrm2(void) const;
    double asum(void) const;
    double amax(void) const;
    double trace(void) const;
    void symmetrize(char uplo);

    // rank-1 update: *this += alpha * x(kx) * (y(ky))^T
    // where x(kx) is row kx of x, and y(ky) is row ky of y
    void ger(double alpha, const DoubleMatrix& x, int kx,
                           const DoubleMatrix& y, int ky);

    // symmetric rank-1 update
    void syr(char uplo, double alpha, const DoubleMatrix& x,
             int ix, char rowcol);

    // get submatrix A(ia:ia+m,ja:ja+n) of A
    void getsub(const DoubleMatrix& a,int m,int n,int ia,int ja);

    // get submatrix A(ia:ia+m,ja:ja+n) of A and store in
    // this(idest:idest+m,jdest:jdest+n)
    void getsub(const DoubleMatrix& a,int m,int n,int isrc,int jsrc,
                int idest, int jdest);

    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(char transa, char transb,
         double alpha, const DoubleMatrix& a,
         const DoubleMatrix& b, double beta);

    // symmetric_matrix * matrix
    // *this = alpha * A * B + beta * this
    void symm(char side, char uplo,
         double alpha, const DoubleMatrix& a,
         const DoubleMatrix& b, double beta);

    // symmetric rank k update
    // this = beta * this + alpha * A * A^T  (trans=='n')
    // this = beta * this + alpha * A^T * A  (trans=='t')
    void syrk(char uplo, char trans,
         double alpha, const DoubleMatrix& a, double beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(double alpha, const DoubleMatrix& a, double beta);
    void transpose(const DoubleMatrix& a);

    void trmm(char side, char uplo, char trans, char diag,
              double alpha, const DoubleMatrix& a);

    // solve triangular system
    void trsm(char side, char uplo, char trans, char diag,
              double alpha, const DoubleMatrix& a);
    void trtrs(char uplo, char trans, char diag, DoubleMatrix& b) const;

    // Cholesky decomposition of a symmetric matrix
    void potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    void potri(char uplo);

    // Polar decomposition, tolerance ||I-X^T*X||<tol or iter<maxiter
    void polar(double tol, int maxiter);

    // LU decomposition
    void lu(std::valarray<int>& ipiv);

    // compute inverse of a square matrix
    void inverse(void);
    // compute inverse and determinant of a square matrix
    double inverse_det(void);
    // compute determinant of a square matrix in LU form
    double det_from_lu(std::valarray<int> ipiv);
    // compute inverse of a square matrix in LU form
    void inverse_from_lu(std::valarray<int>& ipiv);

    // Inverse of triangular matrix
    void trtri(char uplo,char diag);

    double pocon(char) const;
    void sygst(int, char, const DoubleMatrix&);

    // compute eigenvalues and eigenvectors of symmetric matrix *this
    void syev(char uplo, std::valarray<double>& w, DoubleMatrix& z);

    // compute eigenvalues (only) of symmetric matrix *this
    void syev(char uplo, std::valarray<double>& w);

    // compute eigenvalues and eigenvectors of symmetric matrix *this
    // using the divide and conquer method of Tisseur and Dongarra
    void syevd(char uplo, std::valarray<double>& w, DoubleMatrix& z);

    // compute eigenvalues (only) of symmetric matrix *this
    // using the divide and conquer method of Tisseur and Dongarra
    void syevd(char uplo, std::valarray<double>& w);

    // compute eigenvalues and eigenvectors of symmetric matrix *this
    // using the expert driver
    void syevx(char uplo, std::valarray<double>& w, DoubleMatrix& z,
       double abstol);

    // permute the coeff of the matrix *this
    void lapiv(char direc, char rowcol, const int *ipiv);
    // signature of a permutation returned by lu
    int signature(std::valarray<int> ipiv);

    // compute eigenvalues (only) of symmetric matrix *this
    // using the divide and conquer method of Tisseur and Dongarra
    //void syevx(char uplo, std::valarray<double>& w);
};
std::ostream& operator << ( std::ostream& os, const DoubleMatrix& a );

class ComplexMatrix
{
  private:

    const Context& ctxt_;
    int ictxt_;
    int lld_;  // leading dimension of local matrix
    int m_, n_;   // size of global matrix
    int mb_, nb_; // size of blocks
    int mloc_, nloc_;   // size of local array
    int size_; // size of local array;
    int nprow_, npcol_; // number of process rows and cols in  context
    int myrow_, mycol_; // position of my process in the process grid
    int mblocks_, nblocks_; // number of local blocks
    int desc_[9];
    bool active_;
    bool m_incomplete_, n_incomplete_; // this process has an incomplete block
    bool reference_; // object was created using the copy constructor
    std::complex<double>* val;

  public:

    std::complex<double>* valptr(int i=0) { return &val[i]; }
    const std::complex<double>* cvalptr(int i=0) const { return &val[i]; }
    std::complex<double>& operator[] (int i) { return val[i]; }
    const std::complex<double>& operator[] (int i) const { return val[i]; }
    const Context& context(void) const { return ctxt_; }
    int ictxt(void) const { return ictxt_; }
    int lld(void) const { return lld_; }
    int m(void) const { return m_; } // size of global matrix
    int n(void) const { return n_; } // size of global matrix
    int mb(void) const { return mb_; } // size of blocks
    int nb(void) const { return nb_; } // size of blocks
    int mloc(void) const { return mloc_; } // size of local array
    int nloc(void) const { return nloc_; } // size of local array
    int size(void) const { return size_; } // size of local array
    int localsize(void) const { return mloc_*nloc_; } // local size of val
    double memsize(void) const
    { return (double) m_ * (double) n_ * sizeof(std::complex<double>); }
    double localmemsize(void) const
    { return (double) mloc_ * (double) nloc_ * sizeof(std::complex<double>); }
    const int* desc(void) const { return &desc_[0]; }

    // local block size of block (l,m)
    int mbs(int l) const
    {
      // if l is the last block and this process holds an incomplete
      // block, then the size is m_%mb_. Otherwise, it is mb_.
      return ( (l==mblocks_-1) && ( m_incomplete_ ) ) ? m_%mb_ : mb_;
    }
    int nbs(int m) const
    {
      // if m is the last block and this process holds an incomplete
      // block, then the size is n_%nb_. Otherwise, it is nb_.
      return ( (m==nblocks_-1) && ( n_incomplete_ ) ) ? n_%nb_ : nb_;
    }

    // number of local blocks
    int mblocks(void) const { return mblocks_; }
    int nblocks(void) const { return nblocks_; }

    // functions to compute local indices from global indices

    // index of blocks: element (i,j) is at position (x,y)
    // in local block (l,m) of process (pr,pc)
    int l(int i) const { return i/(nprow_*mb_); }
    int x(int i) const { return i % mb_; }
    int pr(int i) const { return (i/mb_) % nprow_; }

    int m(int j) const { return j/(npcol_*nb_); }
    int y(int j) const { return j % nb_; }
    int pc(int j) const { return (j/nb_) % npcol_; }

    // global indices:
    // (i,j) is the global index of element (x,y) of block (l,m)
    int i(int l, int x) const { return (l * nprow_ + myrow_) * mb_ + x; }
    int j(int m, int y) const { return (m * npcol_ + mycol_) * nb_ + y; }

    int iglobal(int ilocal) const
    { return mb_*(nprow_*(ilocal/mb_)+myrow_)+ilocal%mb_; }
    int jglobal(int jlocal) const
    { return nb_*(npcol_*(jlocal/nb_)+mycol_)+jlocal%nb_; }

    int iglobal(int iprow, int ilocal) const
    { return mb_*(nprow_*(ilocal/mb_)+iprow)+ilocal%mb_; }
    int jglobal(int jpcol, int jlocal) const
    { return nb_*(npcol_*(jlocal/nb_)+jpcol)+jlocal%nb_; }

    // store element a(ii,jj) (where ii,jj are global indices)
    // in array val:
    //        int iii = l(ii) * mb_ + x(ii);
    //        int jjj = m(jj) * nb_ + y(jj);
    //        val[iii+mloc_*jjj] = a_ij;

    // active flag: the matrix has elements on this process
    bool active(void) const { return active_; }

    void init_size(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE,
      int nb = MATRIX_DEF_BLOCK_SIZE);

    void resize(int m, int n, int mb = MATRIX_DEF_BLOCK_SIZE,
      int nb = MATRIX_DEF_BLOCK_SIZE)
    {
      const int old_size = size_;
      init_size(m,n,mb,nb);
      if ( size_ == old_size ) return;
      delete[] val;
      val = new std::complex<double>[size_];
      clear();
    }

    void print(std::ostream& os) const;

    explicit ComplexMatrix(const Context& ctxt) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), mloc_(0), nloc_(0),
        size_(0), reference_(false), val(0) {}

    // Construct a ComplexMatrix of dimensions m,n
    explicit ComplexMatrix(const Context& ctxt, int m, int n,
        int mb=MATRIX_DEF_BLOCK_SIZE,
        int nb=MATRIX_DEF_BLOCK_SIZE) : ctxt_(ctxt),
        m_(0), n_(0), mb_(0), nb_(0), size_(0), reference_(false), val(0)
    {
      resize(m,n,mb,nb);
    }

    // copy constructor: create a separate copy of rhs
    explicit ComplexMatrix(const ComplexMatrix& rhs) : ctxt_(rhs.context()),
      reference_(false)
    {
      init_size(rhs.m(),rhs.n(),rhs.mb(),rhs.nb());
      val = new std::complex<double>[size_];
      memcpy(val, rhs.val, size_*sizeof(std::complex<double>));
    }

    // reference constructor: create a proxy for a DoubleMatrix rhs
    explicit ComplexMatrix(DoubleMatrix& rhs);
    explicit ComplexMatrix(const DoubleMatrix& rhs);

    ~ComplexMatrix(void)
    {
      if ( !reference_ ) delete[] val;
    }

    ComplexMatrix& operator=(const ComplexMatrix& a);

    ComplexMatrix& operator+=(const ComplexMatrix& a);
    ComplexMatrix& operator-=(const ComplexMatrix& a);

    ComplexMatrix& operator*=(double a);
    ComplexMatrix& operator*=(std::complex<double> a);

    void matgather(std::complex<double> *a, int lda) const;

    void initdiag(const std::complex<double>* const a); // initialize diagonal
    void init(const std::complex<double>* const a, int lda);
    void set(char uplo, std::complex<double> x);
    void identity(void);
    void clear(void);

    std::complex<double> dot(const ComplexMatrix &a) const;  // tr A^H * A
    std::complex<double> dotu(const ComplexMatrix &a) const; // tr A^T * A
    void axpy(std::complex<double> alpha, const ComplexMatrix &a);
    void axpy(double alpha, const ComplexMatrix &a);
    void scal(std::complex<double> alpha);
    void scal(double alpha);

    double nrm2(void) const;
    double asum(void) const;
    double amax(void) const;
    std::complex<double> trace(void) const;
    void symmetrize(char uplo);

    // rank-1 update: *this += alpha * x(kx) * (y(ky))^T
    // where x(kx) is row kx of x, and y(ky) is row ky of y
    void ger(std::complex<double> alpha, const ComplexMatrix& x,int kx,
                                    const ComplexMatrix& y,int ky);
    void geru(std::complex<double> alpha, const ComplexMatrix& x,int kx,
                                     const ComplexMatrix& y,int ky);
    void gerc(std::complex<double> alpha, const ComplexMatrix& x,int kx,
                                     const ComplexMatrix& y,int ky);

    // symmetric rank-1 update
    void her(char uplo, std::complex<double> alpha,
             const ComplexMatrix& x, int ix, char rowcol);

    // get submatrix A(ia:ia+m,ja:ja+n) of A
    void getsub(const ComplexMatrix& a, int m, int n, int ia, int ja);

    // get submatrix A(ia:ia+m,ja:ja+n) of A and store in
    // this(idest:idest+m,jdest:jdest+n)
    void getsub(const ComplexMatrix& a,int m,int n,int isrc,int jsrc,
                int idest, int jdest);

    // matrix * matrix
    // this = alpha*op(A)*op(B)+beta*this
    void gemm(char transa, char transb,
         std::complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, std::complex<double> beta);

    // hermitian_matrix * matrix
    // *this = alpha * A * B + beta * this
    void hemm(char side, char uplo,
         std::complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, std::complex<double> beta);

    // complex_symmetric_matrix * matrix
    // *this = alpha * A * B + beta * this
    void symm(char side, char uplo,
         std::complex<double> alpha, const ComplexMatrix& a,
         const ComplexMatrix& b, std::complex<double> beta);

    // hermitian rank k update
    // this = beta * this + alpha * A * A^H  (trans=='n')
    // this = beta * this + alpha * A^H * A  (trans=='c')
    void herk(char uplo, char trans,
         double alpha, const ComplexMatrix& a, double beta);

    // matrix transpose
    // this = alpha * transpose(A) + beta * this
    void transpose(std::complex<double> alpha, const ComplexMatrix& a,
                   std::complex<double> beta);
    void transpose(const ComplexMatrix& a);

    void trmm(char side, char uplo, char trans, char diag,
              std::complex<double> alpha, const ComplexMatrix& a);
    void trsm(char side, char uplo, char trans, char diag,
              std::complex<double> alpha, const ComplexMatrix& a);
    void trtrs(char uplo, char trans, char diag, ComplexMatrix& b) const;

    // Cholesky decomposition of a hermitian matrix
    void potrf(char uplo);
    // Inverse of a symmetric matrix from Cholesky factor
    void potri(char uplo);

    // Polar decomposition, tolerance ||I-X^H*X||<tol or iter<maxiter
    void polar(double tol, int maxiter);

    // LU decomposition
    void lu(std::valarray<int>& ipiv);

    // compute inverse of a square matrix
    void inverse(void);
    // compute inverse and determinant of a square matrix
    std::complex<double> inverse_det(void);
    // compute determinant of a square matrix in LU form
    std::complex<double> det_from_lu(std::valarray<int> ipiv);
    // compute inverse of a square matrix in LU form
    void inverse_from_lu(std::valarray<int>& ipiv);

    // Inverse of triangular matrix
    void trtri(char uplo,char diag);

    double pocon(char) const;

    // compute eigenvalues and eigenvectors of hermitian matrix *this
    void heev(char uplo, std::valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues (only) of hermitian matrix *this
    void heev(char uplo, std::valarray<double>& w);

    // eigenvalues/eigenvectors using the divide and conquer method
    // compute eigenvalues and eigenvectors of hermitian matrix *this
    void heevd(char uplo, std::valarray<double>& w, ComplexMatrix& z);
    // compute eigenvalues (only) of hermitian matrix *this
    void heevd(char uplo, std::valarray<double>& w);

    // permute the coeff of the matrix *this
    void lapiv(char direc, char rowcol, const int *ipiv);
    // signature of a permutation returned by lu
    int signature(std::valarray<int> ipiv);
};
std::ostream& operator << ( std::ostream& os, const ComplexMatrix& a );
#endif
