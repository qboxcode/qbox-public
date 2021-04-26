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
// Matrix.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <vector>
#include <complex>
#include <limits>
#include <iostream>
using namespace std;
#ifdef USE_MPI
#include <mpi.h>
#endif
#include "MPIdata.h"

#include "Context.h"
#ifdef SCALAPACK
#include "blacs.h"
#endif

#include "Matrix.h"

#ifdef ADD_
#define numroc     numroc_
#define pdtran     pdtran_
#define pztranc    pztranc_
#define pdsymm     pdsymm_
#define pzsymm     pzsymm_
#define pzhemm     pzhemm_
#define pdgemm     pdgemm_
#define pzgemm     pzgemm_
#define pdsyrk     pdsyrk_
#define pzherk     pzherk_
#define pdsyr      pdsyr_
#define pdger      pdger_
#define pzgerc     pzgerc_
#define pzgeru     pzgeru_
#define pigemr2d   pigemr2d_
#define pdgemr2d   pdgemr2d_
#define pzgemr2d   pzgemr2d_
#define pdtrmm     pdtrmm_
#define pdtrsm     pdtrsm_
#define pztrsm     pztrsm_
#define pdtrtrs    pdtrtrs_
#define pztrtrs    pztrtrs_
#define pdpotrf    pdpotrf_
#define pzpotrf    pzpotrf_
#define pdpotri    pdpotri_
#define pdpocon    pdpocon_
#define pdsygst    pdsygst_
#define pdsyev     pdsyev_
#define pdsyevd    pdsyevd_
#define pdsyevx    pdsyevx_
#define pzheev     pzheev_
#define pzheevd    pzheevd_
#define pdtrtri    pdtrtri_
#define pztrtri    pztrtri_
#define pdlatra    pdlatra_
#define pdlacp2    pdlacp2_
#define pdlacp3    pdlacp3_
#define pdgetrf    pdgetrf_
#define pzgetrf    pzgetrf_
#define pdgetri    pdgetri_
#define pzgetri    pzgetri_
#define pdlapiv    pdlapiv_
#define pzlapiv    pzlapiv_
#define pdlapv2    pdlapv2_
#define pzlapv2    pzlapv2_

#define dscal      dscal_
#define zscal      zscal_
#define zdscal     zdscal_
#define dcopy      dcopy_
#define ddot       ddot_
#define dnrm2      dnrm2_
#define dznrm2     dznrm2_
#define zdotu      zdotu_
#define zdotc      zdotc_
#define daxpy      daxpy_
#define zaxpy      zaxpy_
#define dsymm      dsymm_
#define zsymm      zsymm_
#define zhemm      zhemm_
#define dgemm      dgemm_
#define zgemm      zgemm_
#define dsyr       dsyr_
#define dger       dger_
#define zgerc      zgerc_
#define zgeru      zgeru_
#define dsyrk      dsyrk_
#define zherk      zherk_
#define dtrmm      dtrmm_
#define dtrsm      dtrsm_
#define dtrtri     dtrtri_
#define ztrtri     ztrtri_
#define ztrsm      ztrsm_
#define dtrtrs     dtrtrs_
#define ztrtrs     ztrtrs_
#define dpotrf     dpotrf_
#define zpotrf     zpotrf_
#define dpotri     dpotri_
#define dpocon     dpocon_
#define dsygst     dsygst_
#define dsyev      dsyev_
#define zheev      zheev_
#define idamax     idamax_
#define dgetrf     dgetrf_
#define zgetrf     zgetrf_
#define dgetri     dgetri_
#define zgetri     zgetri_
#endif

extern "C"
{
  int numroc(const int*, const int*, const int*, const int*, const int*);
#ifdef SCALAPACK
  // PBLAS
  void pdsymm(const char*, const char*, const int*, const int*, const double*,
       const double*, const int*, const int*, const int*,
       const double*, const int*, const int*, const int*,
       const double*, double*, const int*, const int*, const int*);
  void pzsymm(const char*, const char*, const int*, const int*,
       const complex<double>*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, complex<double>*, const int*, const int*,
       const int*);
  void pzhemm(const char*, const char*, const int*, const int*,
       const complex<double>*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, complex<double>*, const int*, const int*,
       const int*);
  void pdgemm(const char*, const char*, const int*,
       const int*, const int*, const double*,
       const double*, const int*, const int*, const int*,
       const double*, const int*, const int*, const int*,
       const double*, double*, const int*, const int*, const int*);
  void pzgemm(const char*, const char*, const int*,
       const int*, const int*, const complex<double>*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, const int*, const int*, const int*,
       const complex<double>*, complex<double>*, const int*, const int*,
       const int*);
  void pdger(const int*, const int*, const double*,
       const double*, const int*, const int*, const int*, const int*,
       const double*, const int*, const int*, const int*, const int*,
       double*, const int*, const int*, const int*);
  void pzgerc(const int*, const int*, const complex<double>*,
       const complex<double>*, const int*, const int*, const int*, const int*,
       const complex<double>*, const int*, const int*, const int*, const int*,
       complex<double>*, const int*, const int*, const int*);
  void pzgeru(const int*, const int*, const complex<double>*,
       const complex<double>*, const int*, const int*, const int*, const int*,
       const complex<double>*, const int*, const int*, const int*, const int*,
       complex<double>*, const int*, const int*, const int*);
  void pdsyr(const char*, const int*,
       const double*, const double*, const int*, const int*, const int*,
       const int*, double*, const int*, const int*, const int*);
  void pdsyrk(const char*, const char*, const int*, const int*,
       const double*, const double*, const int*, const int*, const int*,
       const double*, double*, const int*, const int*, const int*);
  void pzherk(const char*, const char*, const int*, const int*,
       const double*, const complex<double>*, const int*,
       const int*, const int*,
       const double*, complex<double>*, const int*,
       const int*, const int*);
  void pdtran(const int*,const  int*, const double*,
       const double*, const int*, const int*, const int*,
       double*, const double*, const int*, const int*, const int*);
  void pztranc(const int*, const int*, const complex<double>*,
       const complex<double>*, const int*, const int*, const int*,
       complex<double>*, const complex<double>*, const int*, const int*,
       const int*);
  void pdtrmm(const char*, const char*, const char*, const char*,
       const int*, const int*, const double*,
       const double*, const int*, const int*, const int*,
       double*, const int*, const int*, const int*);
  void pdtrsm(const char*, const char*, const char*, const char*,
       const int*, const int*, const double*,
       const double*, const int*, const int*, const int*,
       double*, const int*, const int*, const int*);
  void pztrsm(const char*, const char*, const char*, const char*,
       const int*, const int*, const complex<double>*,
       const complex<double>*, const int*, const int*, const int*,
       complex<double>*, const int*, const int*, const int*);
  double pdlatra(const int*,const double*,const int*,const int*,const int*);
  // SCALAPACK
  void pdtrtrs(const char*, const char*, const char*, const int*, const int*,
               const double*, const int*, const int*, const int*,
               double*, const int*, const int*, const int*, int*);
  void pztrtrs(const char*, const char*, const char*, const int*, const int*,
               const complex<double>*, const int*, const int*, const int*,
               complex<double>*, const int*, const int*, const int*, int*);
  void pigemr2d(const int*,const int*,
                const int*,const int*,const int*, const int*,
                int*,const int*,const int*,const int*,const int*);
  void pdgemr2d(const int*,const int*,
                const double*,const int*,const int*, const int*,
                double*,const int*,const int*,const int*,const int*);
  void pzgemr2d(const int*,const int*,
                const complex<double>*,const int*,const int*, const int*,
                complex<double>*,const int*,const int*,const int*,const int*);
  void pdpotrf(const char*, const int*, double*, const int*,
               const int*, const int*, const int*);
  void pzpotrf(const char*, const int*, complex<double>*, const int*,
               const int*, const int*, const int*);
  void pdpotri(const char*, const int*, double*, const int*,
               const int*, const int*, const int*);
  void pdpocon(const char*, const int*, const double*,
               const int*, const int*, const int*, const double*, double*,
               double*, const int*, int*, const int*, int*);
  void pdsygst(const int*, const char*, const int*, double*,
               const int*, const int*, const int*, const double*, const int*,
               const int*, const int*, double*, int*);
  void pdsyev(const char*, const char*, const int*,
              double*, const int*, const int*, const int*, double*, double*,
              const int*, const int*, const int*, double*, const int*, int*);
  void pdsyevd(const char*, const char*, const int*,
              double*, const int*, const int*, const int*, double*, double*,
              const int*, const int*, const int*, double*, const int*, int*,
              int*, int*);
  void pdsyevx(const char* jobz, const char* range, const char* uplo,
               const int* n, double* a, const int* ia, const int* ja,
               const int* desca, double* vl, double* vu,
               const int* il, const int* iu, double* abstol,
               int* nfound, int* nz, double* w,
               const double* orfac, double* z, const int* iz, const int* jz,
               const int* descz, double* work, const int* lwork,
               int* iwork, int* liwork, int* ifail,
               int* icluster, double* gap, int* info);
  void pzheev(const char* jobz, const char* uplo, const int* n,
              complex<double>* a, const int* ia, const int* ja,
              const int* desca, double* w, complex<double> *z,
              const int* iz, const int* jz, const int* descz,
              complex<double>* work, int* lwork,
              double* rwork, int* lrwork, int* info);
  void pzheevd(const char* jobz, const char* uplo, const int* n,
               complex<double>* a, const int* ia, const int* ja,
               const int* desca, double* w, complex<double> *z,
               const int* iz, const int* jz, const int* descz,
               complex<double>* work, int* lwork,
               double* rwork, const int* lrwork,
               int* iwork, int* liwork, int* info);
  void pdtrtri(const char*, const char*, const int*, double*,
               const int*, const int*, const int*, int*);
  void pztrtri(const char*, const char*, const int*, complex<double>*,
               const int*, const int*, const int*, int*);
  void pdgetrf(const int* m, const int* n, double* val,
               int* ia, const int* ja, const int* desca, int* ipiv, int* info);
  void pzgetrf(const int* m, const int* n, complex<double>* val,
               int* ia, const int* ja, const int* desca, int* ipiv, int* info);
  void pdgetri(const int* n, double* val,
               const int* ia, const int* ja, int* desca, int* ipiv,
               double* work, int* lwork, int* iwork, int* liwork, int* info);
  void pzgetri(const int* n, complex<double>* val, const int* ia,
               const int* ja, int* desca, int* ipiv, complex<double>* work,
               int* lwork, int* iwork, int* liwork, int* info);

  void pdlapiv(const char* direc, const char* rowcol, const char* pivroc,
               const int* m, const int* n, double *a, const int* ia,
               const int* ja, const int* desca, int* ipiv, const int* ip,
               const int* jp, const int* descp, int* iwork);
  void pzlapiv(const char* direc, const char* rowcol, const char* pivroc,
               const int* m, const int* n, complex<double> *a, const int* ia,
               const int* ja, const int* desca, int* ipiv, const int* ip,
               const int* jp, const int* descp, int* iwork);
  void pdlapv2(const char* direc, const char *rowcol,
               const int* m, const int *n, double *val,
               const int *ia, const int *ja, const int* desca,
               int *ipiv, const int *ip, const int *jp, const int *descp);
  void pzlapv2(const char* direc, const char *rowcol,
               const int* m, const int *n, complex<double> *val,
               const int *ia, const int *ja, const int* desca,
               int *ipiv, const int *ip, const int *jp, const int *descp);

#endif
  // BLAS1
  void dscal(const int*, const double*, double*, const int*);
  void zscal(const int*, const complex<double>*, complex<double>*, const int*);
  void zdscal(const int*, const double*, complex<double>*, const int*);
  void daxpy(const int *, const double *, const double *, const int *,
             double *, const int *);
  void zaxpy(const int *, const complex<double> *, const complex<double> *,
             const int *, complex<double> *, const int *);
  void dcopy(const int *, const double*, const int *, double*, const int*);
  double ddot(const int *, const double *, const int *,
              const double *, const int *);
  double dnrm2(const int *, const double *, const int *);
  double dznrm2(const int *, const complex<double> *, const int *);
  complex<double> zdotc(const int *, const complex<double>*, const int *,
                        const complex<double>*, const int *);
  complex<double> zdotu(const int *, const complex<double>*, const int *,
                        const complex<double>*, const int *);
  int idamax(const int *, const double*, const int*);
  // BLAS3
  void dsymm(const char*, const char*, const int*, const int *,
             const double*, const double*, const int*,
             const double*, const int*,
             const double*, double*, const int*);
  void zsymm(const char*, const char*, const int*, const int *,
             const complex<double>*, const complex<double>*, const int*,
             const complex<double>*, const int*,
             const complex<double>*, complex<double>*, const int*);
  void zhemm(const char*, const char*, const int*, const int *,
             const complex<double>*, const complex<double>*, const int*,
             const complex<double>*, const int*,
             const complex<double>*, complex<double>*, const int*);
  void dgemm(const char*, const char*, const int*, const int *, const int*,
             const double*, const double*, const int*,
             const double*, const int*,
             const double*, double*, const int*);
  void zgemm(const char*, const char*, const int*, const int *, const int*,
             const complex<double>*, const complex<double>*, const int*,
             const complex<double>*, const int*,
             const complex<double>*, complex<double>*, const int*);
  void zgerc(const int*, const int *, const complex<double>*,
             const complex<double>*, const int*,
             const complex<double>*, const int*,
             const complex<double>*, const int*);
  void zgeru(const int*, const int *, const complex<double>*,
             const complex<double>*, const int*,
             const complex<double>*, const int*,
             const complex<double>*, const int*);
  void dger(const int *, const int*, const double *,
            const double *, const int *, const double *, const int *,
            double*, const int*);
  void dsyr(const char*, const int *, const double *,
            const double *, const int *, double *, const int *);
  void dsyrk(const char*, const char*, const int *, const int *,
             const double *, const double *, const int *,
             const double *, double *, const int *);
  void zherk(const char* uplo, const char* trans, const int* n, const int* k,
             const double* alpha, const complex<double>* a,
             const int*  lda,
             const double* beta, complex<double>* c, const int* ldc);
  void dtrmm(const char*, const char*, const char*, const char*,
             const int*, const int *, const double*, const double*,
             const int*, double*, const int*);
  void dtrsm(const char*, const char*, const char*, const char*,
             const int*, const int *, const double*, const double*,
             const int*, double*, const int*);
  void ztrsm(const char*, const char*, const char*, const char*,
             const int*, const int *, const complex<double>*,
             const complex<double>*, const int*, complex<double>*, const int*);
  // LAPACK
  void dtrtrs(const char*, const char*, const char*,
              const int*, const int*, const double*, const int*,
              double*, const int*, int*);
  void dpotrf(const char*, const int*, double*, const int*, int*);
  void zpotrf(const char*, const int*, complex<double>*, const int*, int*);
  void dpotri(const char*, const int*, double*, const int*, int*);
  void dpocon(const char*, const int *, const double *, const int *,
              const double *, double *, double *, const int *, int *);
  void dsygst(const int*, const char*, const int*,
              double*, const int*, const double*, const int*, int*);
  void dsyev(const char* jobz, const char* uplo, const int* n, double* a,
             const int *lda, double *w, double*work,
             int *lwork, int *info);
  void zheev(const char* jobz, const char* uplo, const int *n,
             complex<double>* a, const int *lda, double* w,
             complex<double>* work, int *lwork, double* rwork, int *info);
  void dtrtri(const char*, const char*, const int*, double*, const int*, int* );
  void ztrtri(const char*, const char*, const int*, complex<double>*,
              const int*, int* );
  void ztrtrs(const char*, const char*, char*, const int*, const int*,
              complex<double>*, const int*, complex<double>*, int*, int* );
  void dgetrf(const int* m, const int* n, double* a, const int* lda,
              int* ipiv, int*info);
  void zgetrf(const int* m, const int* n, complex<double>* a, const int* lda,
              int* ipiv, int*info);
  void dgetri(const int* m, double* val, const int* lda, int* ipiv,
              double* work, int* lwork, int* info);
  void zgetri(const int* m, complex<double>* val, const int* lda, int* ipiv,
              complex<double>* work, int* lwork, int* info);
}

////////////////////////////////////////////////////////////////////////////////
// numroc0: ScaLAPACK numroc function specialized for the case isrcproc=0
// i.e. the process holding the first row/col of the matrix is proc 0
int numroc0(int n, int nb, int iproc, int nprocs)
{
  // n       number of rows/cols of the distributed matrix
  // nb      block size
  // iproc   coordinate of the process whose local array size is being computed
  //         iproc = 0..nprocs
  // nprocs  number of processes over which the matrix is distributed

  // nblocks = total number of whole nb blocks
  int n_whole_blocks = n / nb;

  // minimum number of rows or cols a process can have
  int nroc = ( n_whole_blocks / nprocs ) * nb;

  // number of extra blocks needed
  int n_extra_blocks = n_whole_blocks % nprocs;

  // adjust numroc depending on iproc
  if ( iproc < n_extra_blocks )
    nroc += nb;
  else if ( iproc == n_extra_blocks )
    nroc += n % nb;

  return nroc;
}

////////////////////////////////////////////////////////////////////////////////
int DoubleMatrix::mloc(int irow) const
{
  return numroc0(m_,mb_,irow,nprow_);
}

////////////////////////////////////////////////////////////////////////////////
int DoubleMatrix::nloc(int icol) const
{
  return numroc0(n_,nb_,icol,npcol_);
}

////////////////////////////////////////////////////////////////////////////////
// reference constructor create a proxy for a ComplexMatrix rhs
DoubleMatrix::DoubleMatrix(ComplexMatrix& rhs) : ctxt_(rhs.context()),
  reference_(true)
{
  int new_m = 2 * rhs.m();
  int new_mb = 2 * rhs.mb();
  init_size(new_m,rhs.n(),new_mb,rhs.nb());
  val = (double*) rhs.valptr();
}

////////////////////////////////////////////////////////////////////////////////
// reference constructor create a proxy for a const ComplexMatrix rhs
DoubleMatrix::DoubleMatrix(const ComplexMatrix& rhs) : ctxt_(rhs.context()),
  reference_(true)
{
  int new_m = 2 * rhs.m();
  int new_mb = 2 * rhs.mb();
  init_size(new_m,rhs.n(),new_mb,rhs.nb());
  val = (double*) rhs.cvalptr();
}

////////////////////////////////////////////////////////////////////////////////
// reference constructor create a proxy for a DoubleMatrix rhs
ComplexMatrix::ComplexMatrix(DoubleMatrix& rhs) : ctxt_(rhs.context()),
  reference_(true)
{
  assert(rhs.m()%2 == 0);
  int new_m = rhs.m() / 2;
  assert(rhs.mb()%2 == 0);
  int new_mb = rhs.mb() / 2;
  init_size(new_m,rhs.n(),new_mb,rhs.nb());
  val = (complex<double>*) rhs.valptr();
}

////////////////////////////////////////////////////////////////////////////////
// reference constructor create a proxy for a const DoubleMatrix rhs
ComplexMatrix::ComplexMatrix(const DoubleMatrix& rhs) : ctxt_(rhs.context()),
  reference_(true)
{
  assert(rhs.m()%2 == 0);
  int new_m = rhs.m() / 2;
  assert(rhs.mb()%2 == 0);
  int new_mb = rhs.mb() / 2;
  init_size(new_m,rhs.n(),new_mb,rhs.nb());
  val = (complex<double>*) rhs.cvalptr();
}

////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::init_size(int m, int n, int mb, int nb)
{
  assert(m>=0);
  assert(n>=0);
  assert(mb>=0);
  assert(nb>=0);
  m_ = m;
  n_ = n;
#ifdef SCALAPACK
  mb_ = mb;
  nb_ = nb;
#else
  mb_ = m;
  nb_ = n;
#endif
  if ( mb_ == 0 ) mb_ = 1;
  if ( nb_ == 0 ) nb_ = 1;
  ictxt_ = ctxt_.ictxt();
  nprow_ = ctxt_.nprow();
  npcol_ = ctxt_.npcol();
  myrow_ = ctxt_.myrow();
  mycol_ = ctxt_.mycol();
  active_ = myrow_ >= 0;
  int isrcproc=0;
  mloc_ = 0;
  nloc_ = 0;
  if ( m_ != 0 )
    mloc_ = numroc(&m_,&mb_,&myrow_,&isrcproc,&nprow_);
  if ( n_ != 0 )
    nloc_ = numroc(&n_,&nb_,&mycol_,&isrcproc,&npcol_);
  size_ = mloc_ * nloc_;

  // set leading dimension of val array to mloc_;
  lld_ = mloc_;
  if ( lld_ == 0 ) lld_ = 1;

  // total and local number of blocks
  mblocks_ = 0;
  nblocks_ = 0;
  m_incomplete_ = false;
  n_incomplete_ = false;
  if ( active_ && mb_ > 0 && nb_ > 0 )
  {
    mblocks_ = ( mloc_ + mb_ - 1 ) / mb_;
    nblocks_ = ( nloc_ + nb_ - 1 ) / nb_;
    m_incomplete_ = mloc_ % mb_ != 0;
    n_incomplete_ = nloc_ % nb_ != 0;
  }

  if ( active_ )
  {
    desc_[0] = 1;
  }
  else
  {
    desc_[0] = -1;
  }
  desc_[1] = ictxt_;
  desc_[2] = m_;
  desc_[3] = n_;
  desc_[4] = mb_;
  desc_[5] = nb_;
  desc_[6] = 0;
  desc_[7] = 0;
  desc_[8] = lld_;
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::init_size(int m, int n, int mb, int nb)
{
  assert(m>=0);
  assert(n>=0);
  assert(mb>=0);
  assert(nb>=0);
  m_ = m;
  n_ = n;
#ifdef SCALAPACK
  mb_ = mb;
  nb_ = nb;
#else
  mb_ = m;
  nb_ = n;
#endif
  if ( mb_ == 0 ) mb_ = 1;
  if ( nb_ == 0 ) nb_ = 1;
  ictxt_ = ctxt_.ictxt();
  nprow_ = ctxt_.nprow();
  npcol_ = ctxt_.npcol();
  myrow_ = ctxt_.myrow();
  mycol_ = ctxt_.mycol();
  active_ = myrow_ >= 0;
  int isrcproc=0;
  mloc_ = 0;
  nloc_ = 0;

  if ( m_ != 0 )
    mloc_ = numroc(&m_,&mb_,&myrow_,&isrcproc,&nprow_);
  if ( n_ != 0 )
    nloc_ = numroc(&n_,&nb_,&mycol_,&isrcproc,&npcol_);
  size_ = mloc_ * nloc_;

  // set leading dimension of val array to mloc_;
  lld_ = mloc_;
  if ( lld_ == 0 ) lld_ = 1;

  // total and local number of blocks
  mblocks_ = 0;
  nblocks_ = 0;
  m_incomplete_ = false;
  n_incomplete_ = false;
  if ( active_ && mb_ > 0 && nb_ > 0 )
  {
    mblocks_ = ( mloc_ + mb_ - 1 ) / mb_;
    nblocks_ = ( nloc_ + nb_ - 1 ) / nb_;
    m_incomplete_ = mloc_ % mb_ != 0;
    n_incomplete_ = nloc_ % nb_ != 0;
  }

  if ( active_ )
  {
    desc_[0] = 1;
  }
  else
  {
    desc_[0] = -1;
  }
  desc_[1] = ictxt_;
  desc_[2] = m_;
  desc_[3] = n_;
  desc_[4] = mb_;
  desc_[5] = nb_;
  desc_[6] = 0;
  desc_[7] = 0;
  desc_[8] = lld_;
}

////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::clear(void)
{
  assert(val!=0||size_==0);
  memset(val,0,size_*sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::clear(void)
{
  assert(val!=0||size_==0);
  memset(val,0,size_*sizeof(complex<double>));
}

////////////////////////////////////////////////////////////////////////////////
// real identity: initialize matrix to identity
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::identity(void)
{
  clear();
  set('d',1.0);
}

////////////////////////////////////////////////////////////////////////////////
// complex identity: initialize matrix to identity
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::identity(void)
{
  clear();
  set('d',complex<double>(1.0,0.0));
}

////////////////////////////////////////////////////////////////////////////////
// set value of diagonal or off-diagonal elements to a constant
// uplo=='u': set strictly upper part to x
// uplo=='l': set strictly lower part to x
// uplo=='d': set diagonal to x
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::set(char uplo, double xx)
{
  if ( active_ )
  {
    if ( uplo=='l' || uplo=='L' )
    {
      // initialize strictly lower part
      for (int li=0; li < mblocks_;li++)
      {
        for (int lj=0; lj < nblocks_;lj++)
        {
          for (int ii=0; ii < mbs(li); ii++)
          {
            for (int jj=0; jj < nbs(lj);jj++)
            {
              if ( i(li,ii) > j(lj,jj) )
                val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ] = xx;
            }
          }
        }
      }
    }
    else if ( uplo=='u' || uplo=='U' )
    {
      // initialize strictly upper part
      for ( int li=0; li < mblocks_; li++ )
      {
        for ( int lj=0; lj < nblocks_; lj++ )
        {
          for ( int ii=0; ii < mbs(li); ii++ )
          {
            for ( int jj=0; jj < nbs(lj); jj++ )
            {
              if ( i(li,ii) < j(lj,jj) )
                val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ] = xx;
            }
          }
        }
      }
    }
    else if ( uplo=='d' || uplo=='D' )
    {
      // initialize diagonal elements
      if ( active() )
      {
        // loop through all local blocks (ll,mm)
        for ( int ll = 0; ll < mblocks(); ll++)
        {
          for ( int mm = 0; mm < nblocks(); mm++)
          {
            // check if block (ll,mm) has diagonal elements
            int imin = i(ll,0);
            int imax = imin + mbs(ll)-1;
            int jmin = j(mm,0);
            int jmax = jmin + nbs(mm)-1;
            // cout << " process (" << myrow_ << "," << mycol_ << ")"
            // << " block (" << ll << "," << mm << ")"
            // << " imin/imax=" << imin << "/" << imax
            // << " jmin/jmax=" << jmin << "/" << jmax << endl;

            if ((imin <= jmax) && (imax >= jmin))
            {
              // block (ll,mm) holds diagonal elements
              int idiagmin = max(imin,jmin);
              int idiagmax = min(imax,jmax);

              // cout << " process (" << myrow_ << "," << mycol_ << ")"
              // << " holds diagonal elements " << idiagmin << " to " <<
              // idiagmax << " in block (" << ll << "," << mm << ")" << endl;

              for ( int ii = idiagmin; ii <= idiagmax; ii++ )
              {
                // access element (ii,ii)
                int jj = ii;
                int iii = ll * mb_ + x(ii);
                int jjj = mm * nb_ + y(jj);
                val[iii+mloc_*jjj] = xx;
              }
            }
          }
        }
      }
    }
    else
    {
      cout << " DoubleMatrix::set: invalid argument" << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::set(char uplo, complex<double> xx)
{
  if ( active_ )
  {
    if ( uplo=='l' || uplo=='L' )
    {
      // initialize strictly lower part
      for (int li=0; li < mblocks_;li++)
      {
        for (int lj=0; lj < nblocks_;lj++)
        {
          for (int ii=0; ii < mbs(li); ii++)
          {
            for (int jj=0; jj < nbs(lj);jj++)
            {
              if ( i(li,ii) > j(lj,jj) )
                val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ] = xx;
            }
          }
        }
      }
    }
    else if ( uplo=='u' || uplo=='U' )
    {
      // initialize strictly upper part
      for ( int li=0; li < mblocks_; li++ )
      {
        for ( int lj=0; lj < nblocks_; lj++ )
        {
          for ( int ii=0; ii < mbs(li); ii++ )
          {
            for ( int jj=0; jj < nbs(lj); jj++ )
            {
              if ( i(li,ii) < j(lj,jj) )
                val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ] = xx;
            }
          }
        }
      }
    }
    else if ( uplo=='d' || uplo=='D' )
    {
      // initialize diagonal elements
      if ( active() )
      {
        // loop through all local blocks (ll,mm)
        for ( int ll = 0; ll < mblocks(); ll++)
        {
          for ( int mm = 0; mm < nblocks(); mm++)
          {
            // check if block (ll,mm) has diagonal elements
            int imin = i(ll,0);
            int imax = imin + mbs(ll)-1;
            int jmin = j(mm,0);
            int jmax = jmin + nbs(mm)-1;
            // cout << " process (" << myrow_ << "," << mycol_ << ")"
            // << " block (" << ll << "," << mm << ")"
            // << " imin/imax=" << imin << "/" << imax
            // << " jmin/jmax=" << jmin << "/" << jmax << endl;

            if ((imin <= jmax) && (imax >= jmin))
            {
              // block (ll,mm) holds diagonal elements
              int idiagmin = max(imin,jmin);
              int idiagmax = min(imax,jmax);

              // cout << " process (" << myrow_ << "," << mycol_ << ")"
              // << " holds diagonal elements " << idiagmin << " to " <<
              // idiagmax << " in block (" << ll << "," << mm << ")" << endl;

              for ( int ii = idiagmin; ii <= idiagmax; ii++ )
              {
                // access element (ii,ii)
                int jj = ii;
                int iii = ll * mb_ + x(ii);
                int jjj = mm * nb_ + y(jj);
                val[iii+mloc_*jjj] = xx;
              }
            }
          }
        }
      }
    }
    else
    {
      cout << " DoubleMatrix::set: invalid argument" << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// initialize *this using a replicated matrix a
void DoubleMatrix::init(const double* const a, int lda)
{
  if ( active_ )
  {
    for ( int li=0; li < mblocks_; li++ )
    {
      for ( int lj=0; lj < nblocks_; lj++ )
      {
        for ( int ii=0; ii < mbs(li); ii++ )
        {
          for ( int jj=0; jj < nbs(lj); jj++ )
          {
            val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ]
                = a[ i(li,ii) + j(lj,jj)*lda ];
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::dot(const DoubleMatrix &x) const
{
  assert( ictxt_ == x.ictxt() );
  double  sum=0.;
  double  tsum=0.;
  if ( active_ )
  {
    assert( m_ == x.m() );
    assert( n_ == x.n() );
    assert( mb_ == x.mb() );
    assert( nb_ == x.nb() );
    assert( mloc_ == x.mloc() );
    assert( nloc_ == x.nloc() );
    assert(size_==x.size());
    int ione=1;
    tsum=ddot(&size_, val, &ione, x.val, &ione);
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
#else
  sum=tsum;
#endif
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
complex<double> ComplexMatrix::dot(const ComplexMatrix &x) const
{
  assert( ictxt_ == x.ictxt() );
  complex<double>  sum=0.0;
  complex<double>  tsum=0.0;
  if ( active_ )
  {
    assert( m_ == x.m() );
    assert( n_ == x.n() );
    assert( mb_ == x.mb() );
    assert( nb_ == x.nb() );
    assert( mloc_ == x.mloc() );
    assert( nloc_ == x.nloc() );
    assert(size_==x.size());
    //int ione=1;
    //tsum=zdotc(&size_, val, &ione, x.val, &ione);
    for ( int i = 0; i < size_; i++ )
      tsum += conj(val[i]) * x.val[i];
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce((double*)&tsum, (double*)&sum, 2,
                  MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
#else
  sum=tsum;
#endif
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
complex<double> ComplexMatrix::dotu(const ComplexMatrix &x) const
{
  assert( ictxt_ == x.ictxt() );
  complex<double>  sum=0.0;
  complex<double>  tsum=0.0;
  if ( active_ )
  {
    assert( m_ == x.m() );
    assert( n_ == x.n() );
    assert( mb_ == x.mb() );
    assert( nb_ == x.nb() );
    assert( mloc_ == x.mloc() );
    assert( nloc_ == x.nloc() );
    assert(size_==x.size());
    //int ione=1;
    //tsum=zdotu(&size_, val, &ione, x.val, &ione);
    for ( int i = 0; i < size_; i++ )
      tsum += val[i] * x.val[i];
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce((double*)&tsum, (double*)&sum, 2,
                  MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
#else
  sum=tsum;
#endif
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::amax(void) const
{
  double am = 0.0, tam = 0.0;
  if ( active_ )
  {
    int ione=1;
    tam = val[idamax(&size_,val,&ione) - 1];
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce(&tam, &am, 1, MPI_DOUBLE, MPI_MAX, ctxt_.comm() );
#else
  am=tam;
#endif
  return am;
}

////////////////////////////////////////////////////////////////////////////////
// axpy: *this = *this + alpha * x
void DoubleMatrix::axpy(double alpha, const DoubleMatrix &x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  if( active_ )
    daxpy(&size_, &alpha, x.val, &ione, val, &ione);
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::axpy(complex<double> alpha, const ComplexMatrix &x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  if( active_ )
    zaxpy(&size_, &alpha, x.val, &ione, val, &ione);
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::axpy(double alpha, const ComplexMatrix &x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  int len = 2 * size_;
  if( active_ )
    daxpy(&len, &alpha, (double*) x.val, &ione, (double*) val, &ione);
}

////////////////////////////////////////////////////////////////////////////////
// real getsub: *this = sub(A)
// copy submatrix A(ia:ia+m, ja:ja+n) into *this;
// *this and A may live in different contexts
void DoubleMatrix::getsub(const DoubleMatrix &a,
  int m, int n, int ia, int ja)
{
#if SCALAPACK
  int iap=ia+1;
  int jap=ja+1;
  assert(n<=n_);
  assert(n<=a.n());
  assert(m<=m_);
  assert(m<=a.m());
  int ione = 1;
  int gictxt;
  Cblacs_get( 0, 0, &gictxt );
  pdgemr2d(&m,&n,a.val,&iap,&jap,a.desc_,val,&ione,&ione,desc_,&gictxt);
#else
  for ( int j = 0; j < n; j++ )
    for ( int i = 0; i < m; i++ )
      val[i+j*m_] = a.val[(i+ia) + (j+ja)*a.m()];
#endif
}

////////////////////////////////////////////////////////////////////////////////
// real getsub: *this = sub(A)
// copy submatrix A(ia:ia+m, ja:ja+n) into *this(idest:idest+m,jdest:jdest+n)
// *this and A may live in different contexts
void DoubleMatrix::getsub(const DoubleMatrix &a,
  int m, int n, int isrc, int jsrc, int idest, int jdest)
{
#if SCALAPACK
  int iap=isrc+1;
  int jap=jsrc+1;
  int idp=idest+1;
  int jdp=jdest+1;
  assert(n<=n_);
  assert(n<=a.n());
  assert(m<=m_);
  assert(m<=a.m());
  int gictxt;
  Cblacs_get( 0, 0, &gictxt );
  pdgemr2d(&m,&n,a.val,&iap,&jap,a.desc_,val,&idp,&jdp,desc_,&gictxt);
#else
  for ( int j = 0; j < n; j++ )
    for ( int i = 0; i < m; i++ )
      val[(idest+i)+(jdest+j)*m_] = a.val[(i+isrc) + (j+jsrc)*a.m()];
#endif
}

////////////////////////////////////////////////////////////////////////////////
// complex getsub: *this = sub(A)
// copy submatrix A(ia:ia+m, ja:ja+n) into *this
// *this and A may live in different contexts
void ComplexMatrix::getsub(const ComplexMatrix &a,
  int m, int n, int ia, int ja)
{
#if SCALAPACK
  int iap=ia+1;
  int jap=ja+1;
  assert(n<=n_);
  assert(n<=a.n());
  assert(m<=m_);
  assert(m<=a.m());
  int ione = 1;
  int gictxt;
  Cblacs_get( 0, 0, &gictxt );
  pzgemr2d(&m,&n,a.val,&iap,&jap,a.desc_,val,&ione,&ione,desc_,&gictxt);
#else
  for ( int j = 0; j < n; j++ )
    for ( int i = 0; i < m; i++ )
      val[i+j*m_] = a.val[(i+ia) + (j+ja)*a.m()];
#endif
}

////////////////////////////////////////////////////////////////////////////////
// complex getsub: *this = sub(A)
// copy submatrix A(ia:ia+m, ja:ja+n) into *this(idest:idest+m,jdest:jdest+n)
// *this and A may live in different contexts
void ComplexMatrix::getsub(const ComplexMatrix &a,
  int m, int n, int isrc, int jsrc, int idest, int jdest)
{
#if SCALAPACK
  int iap=isrc+1;
  int jap=jsrc+1;
  int idp=idest+1;
  int jdp=jdest+1;
  assert(n<=n_);
  assert(n<=a.n());
  assert(m<=m_);
  assert(m<=a.m());
  int gictxt;
  Cblacs_get( 0, 0, &gictxt );
  pzgemr2d(&m,&n,a.val,&iap,&jap,a.desc_,val,&idp,&jdp,desc_,&gictxt);
#else
  for ( int j = 0; j < n; j++ )
    for ( int i = 0; i < m; i++ )
      val[(idest+i)+(jdest+j)*m_] = a.val[(i+isrc) + (j+jsrc)*a.m()];
#endif
}

////////////////////////////////////////////////////////////////////////////////
// real matrix transpose
// this = alpha * transpose(A) + beta * this
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::transpose(double alpha, const DoubleMatrix& a, double beta)
{
  assert(this != &a);
  assert( ictxt_ == a.ictxt() );

  if ( active() )
  {
    assert(a.m() == n_);
    assert(a.n() == m_);

#ifdef SCALAPACK
    int ione = 1;
    pdtran(&m_, &n_, &alpha,
         a.val, &ione, &ione, a.desc_,
         &beta, val, &ione, &ione, desc_);
#else
    scal(beta);
    for ( int i=0; i<m_; i++ )
      for ( int j=0; j<i; j++ )
      {
        val[i*m_+j] += alpha * a.val[j*m_+i];
        val[j*m_+i] += alpha * a.val[i*m_+j];
      }
    for ( int i=0; i<m_; i++ )
      val[i*m_+i] += alpha * a.val[i*m_+i];
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
// real matrix transpose
// *this = transpose(a)
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::transpose(const DoubleMatrix& a)
{
  assert(this != &a);
  transpose(1.0,a,0.0);
}

////////////////////////////////////////////////////////////////////////////////
// complex hermitian transpose
// this = alpha * A^H + beta * this
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::transpose(complex<double> alpha, const ComplexMatrix& a,
  complex<double> beta)
{
  assert(this != &a);
  assert( ictxt_ == a.ictxt() );

  if ( active() )
  {
    assert(a.m() == n_);
    assert(a.n() == m_);

#ifdef SCALAPACK
    int ione = 1;
    pztranc(&m_, &n_, &alpha,
         a.val, &ione, &ione, a.desc_,
         &beta, val, &ione, &ione, desc_);
#else
    scal(beta);
    for ( int i=0; i<m_; i++ )
      for ( int j=0; j<i; j++ )
      {
        val[i*m_+j] += alpha * conj(a.val[j*m_+i]);
        val[j*m_+i] += alpha * conj(a.val[i*m_+j]);
      }
    for ( int i=0; i<m_; i++ )
      val[i*m_+i] += alpha * a.val[i*m_+i];
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
// complex matrix transpose
// *this = transpose(a)
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::transpose(const ComplexMatrix& a)
{
  assert(this != &a);
  transpose(complex<double>(1.0,0.0),a,complex<double>(0.0,0.0));
}

////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::symmetrize(char uplo)
{
  // symmetrize
  // if uplo == 'l' : copy strictly lower triangle to strictly upper triangle
  // if uplo == 'u' : copy strictly upper triangle to strictly lower triangle
  // if uplo == 'n' : A = 0.5 * ( A^T + A )

  if ( uplo == 'n' )
  {
    DoubleMatrix tmp(*this);
    transpose(0.5,tmp,0.5);
  }
  else if ( uplo == 'l' )
  {
    set('u',0.0);
    DoubleMatrix tmp(*this);
    tmp.set('d',0.0);
    transpose(1.0,tmp,1.0);
  }
  else if ( uplo == 'u' )
  {
    set('l',0.0);
    DoubleMatrix tmp(*this);
    tmp.set('d',0.0);
    transpose(1.0,tmp,1.0);
  }
  else
  {
    cout << " DoubleMatrix::symmetrize: invalid argument" << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 2);
#else
      exit(2);
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::symmetrize(char uplo)
{
  // symmetrize
  // uplo == 'l' : copy conjugate of strictly lower triangle to strictly upper
  // uplo == 'u' : copy conjugate of strictly upper triangle to strictly lower
  // uplo == 'n' : A = 0.5 * ( A^H + A )

  if ( uplo == 'n' )
  {
    ComplexMatrix tmp(*this);
    transpose(complex<double>(0.5,0.0),tmp,complex<double>(0.5,0.0));
  }
  else if ( uplo == 'l' )
  {
    set('u',complex<double>(0.0,0.0));
    ComplexMatrix tmp(*this);
    tmp.set('d',complex<double>(0.0,0.0));
    transpose(complex<double>(1.0,0.0),tmp,complex<double>(1.0,0.0));
  }
  else if ( uplo == 'u' )
  {
    set('l',complex<double>(0.0,0.0));
    ComplexMatrix tmp(*this);
    tmp.set('d',complex<double>(0.0,0.0));
    transpose(complex<double>(1.0,0.0),tmp,complex<double>(1.0,0.0));
  }
  else
  {
    cout << " ComplexMatrix::symmetrize: invalid argument" << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 2);
#else
      exit(2);
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::nrm2(void) const
{
  double  sum=0.;
  double  tsum=0.;
  if ( active_ )
  {
    int ione=1;
    // dnrm2 returns sqrt(sum_i x[i]*x[i])
    tsum = dnrm2(&size_,val,&ione);
    tsum = tsum*tsum;
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
#else
  sum=tsum;
#endif
  return sqrt(sum);
}

////////////////////////////////////////////////////////////////////////////////
double ComplexMatrix::nrm2(void) const
{
  double  sum=0.;
  double  tsum=0.;
  if ( active_ )
  {
    int ione=1;
    // dznrm2 returns sqrt(sum_i conjg(x[i])*x[i])
    tsum = dznrm2(&size_,val,&ione);
    tsum = tsum*tsum;
  }
#ifdef SCALAPACK
  if ( active_ )
    MPI_Allreduce(&tsum, &sum, 1, MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
#else
  sum=tsum;
#endif
  return sqrt(sum);
}

////////////////////////////////////////////////////////////////////////////////
// rank-1 update using row kx of x and (row ky of y)^T
// *this = *this + alpha * x(kx) * y(ky)^T
void DoubleMatrix::ger(double alpha,
  const DoubleMatrix& x, int kx, const DoubleMatrix& y, int ky)
{
  assert(x.n()==m_);
  assert(y.n()==n_);
#if SCALAPACK
  int ione=1;

  int ix = kx+1;
  int jx = 1;
  int incx = x.m();

  int iy = ky+1;
  int jy = 1;
  int incy = y.m();
  pdger(&m_,&n_,&alpha,x.val,&ix,&jx,x.desc_,&incx,
                       y.val,&iy,&jy,y.desc_,&incy,
                       val,&ione,&ione,desc_);
#else
  int incx = x.m();
  int incy = y.m();
  dger(&m_,&n_,&alpha,&x.val[kx*x.m()],&incx,
                      &y.val[ky*y.m()],&incy,val,&m_);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// rank-1 update using row kx of x and conj(row ky of y)^T
// *this = *this + alpha * x(kx) * conj(y(ky))^T
void ComplexMatrix::gerc(complex<double> alpha,
  const ComplexMatrix& x, int kx, const ComplexMatrix& y, int ky)
{
  assert(x.n()==m_);
  assert(y.n()==n_);
#if SCALAPACK
  int ione=1;

  int ix = kx+1;
  int jx = 1;
  int incx = x.m();

  int iy = ky+1;
  int jy = 1;
  int incy = y.m();
  pzgerc(&m_,&n_,&alpha,x.val,&ix,&jx,x.desc_,&incx,
                       y.val,&iy,&jy,y.desc_,&incy,
                       val,&ione,&ione,desc_);
#else
  int incx = x.m();
  int incy = y.m();
  zgerc(&m_,&n_,&alpha,&x.val[kx*x.m()],&incx,
                      &y.val[ky*y.m()],&incy,val,&m_);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// rank-1 update using row kx of x and conj(row ky of y)^T
// *this = *this + alpha * x(kx) * y(ky)^T
void ComplexMatrix::geru(complex<double> alpha,
  const ComplexMatrix& x, int kx, const ComplexMatrix& y, int ky)
{
  assert(x.n()==m_);
  assert(y.n()==n_);
#if SCALAPACK
  int ione=1;

  int ix = kx+1;
  int jx = 1;
  int incx = x.m();

  int iy = ky+1;
  int jy = 1;
  int incy = y.m();
  pzgeru(&m_,&n_,&alpha,x.val,&ix,&jx,x.desc_,&incx,
                       y.val,&iy,&jy,y.desc_,&incy,
                       val,&ione,&ione,desc_);
#else
  int incx = x.m();
  int incy = y.m();
  zgeru(&m_,&n_,&alpha,&x.val[kx*x.m()],&incx,
                      &y.val[ky*y.m()],&incy,val,&m_);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// symmetric rank-1 update using a row or a column of a Matrix x
void DoubleMatrix::syr(char uplo, double alpha,
  const DoubleMatrix& x, int k, char rowcol)
{
  assert(n_==m_);
#if SCALAPACK
  int ix,jx,incx,ione=1;
  if ( rowcol == 'c' )
  {
    // use column k of matrix x
    assert(x.m()==n_);
    ix = 1;
    jx = k+1;
    incx = 1;
  }
  else if ( rowcol == 'r' )
  {
    // use row k of matrix x
    assert(x.n()==n_);
    ix = k+1;
    jx = 1;
    incx = x.m();
  }
  else
  {
    cout << " DoubleMatrix::syr: invalid argument rowcol" << endl;
    MPI_Abort(MPI_COMM_WORLD,2);
  }
  pdsyr(&uplo,&n_,&alpha,x.val,&ix,&jx,x.desc_,&incx,
        val,&ione,&ione,desc_);
#else
  if ( rowcol == 'c' )
  {
    // use column k of matrix x
    assert(x.m()==n_);
    int incx = 1;
    dsyr(&uplo,&n_,&alpha,&x.val[k*x.m()],&incx,val,&m_);
  }
  else if ( rowcol == 'r' )
  {
    // use row k of matrix x
    assert(x.n()==n_);
    int incx = x.m();
    dsyr(&uplo,&n_,&alpha,&x.val[k],&incx,val,&m_);
  }
  else
  {
    cout << " DoubleMatrix::syr: invalid argument rowcol" << endl;
    exit(2);
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
DoubleMatrix& DoubleMatrix::operator=(const DoubleMatrix& a)
{
  if ( this == &a ) return *this;

  // operator= works only for matrices having same distribution on same context
  assert( a.ictxt() == ictxt_ && a.m() == m_ && a.mb() == mb_ &&
          a.n() == n_ && a.nb() == nb_ );
  if ( active() )
  {
    for ( int i = 0; i < 9; i++ )
    {
      assert( desc_[i] == a.desc_[i] );
    }
    memcpy(val, a.val, mloc_*nloc_*sizeof(double));
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
ComplexMatrix& ComplexMatrix::operator=(const ComplexMatrix& a)
{
  if ( this == &a ) return *this;

  assert( a.ictxt() == ictxt_ && a.m() == m_ && a.mb() == mb_ &&
          a.n() == n_ && a.nb() == nb_ );
  if ( active() )
  {
    for ( int i = 0; i < 9; i++ )
    {
      assert( desc_[i] == a.desc_[i] );
    }
    memcpy(val, a.val, mloc_*nloc_*sizeof(complex<double>));
  }
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator+=
DoubleMatrix& DoubleMatrix::operator+=(const DoubleMatrix &x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  double alpha = 1.0;
  if( active_ )
    daxpy(&size_, &alpha, x.val, &ione, val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator-=
DoubleMatrix& DoubleMatrix::operator-=(const DoubleMatrix &x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  double alpha = -1.0;
  if( active_ )
    daxpy(&size_, &alpha, x.val, &ione, val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator+=
ComplexMatrix& ComplexMatrix::operator+=(const ComplexMatrix& x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  double alpha = 1.0;
  int two_size = 2 * size_;
  if( active_ )
    daxpy(&two_size, &alpha, (double*) x.val, &ione, (double*) val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator-=
ComplexMatrix& ComplexMatrix::operator-=(const ComplexMatrix& x)
{
  assert( ictxt_ == x.ictxt() );
  int ione=1;
  assert(m_==x.m());
  assert(n_==x.n());
  assert(mloc_==x.mloc());
  assert(nloc_==x.nloc());
  double alpha = -1.0;
  int two_size = 2 * size_;
  if( active_ )
    daxpy(&two_size, &alpha, (double*) x.val, &ione, (double*) val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator*=
DoubleMatrix& DoubleMatrix::operator*=(double alpha)
{
  int ione=1;
  if( active_ )
    dscal(&size_, &alpha, val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
ComplexMatrix& ComplexMatrix::operator*=(double alpha)
{
  int ione=1;
  if( active_ )
    zdscal(&size_, &alpha, val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// operator*=
ComplexMatrix& ComplexMatrix::operator*=(complex<double> alpha)
{
  int ione=1;
  if( active_ )
    zscal(&size_, &alpha, val, &ione);
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
// scal
void DoubleMatrix::scal(double alpha)
{
  *this *= alpha;
}

////////////////////////////////////////////////////////////////////////////////
// scal
void ComplexMatrix::scal(double alpha)
{
  *this *= alpha;
}

////////////////////////////////////////////////////////////////////////////////
// scal
void ComplexMatrix::scal(complex<double> alpha)
{
  *this *= alpha;
}

////////////////////////////////////////////////////////////////////////////////
// matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::gemm(char transa, char transb,
                double alpha, const DoubleMatrix& a,
                const DoubleMatrix& b, double beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( ictxt_ == b.ictxt() );

  if ( active() )
  {
    int m, n, k;
    if ( transa == 'N' || transa == 'n' )
    {
      m = a.m();
      k = a.n();
      assert(a.m()==m_);
    }
    else
    {
      m = a.n();
      k = a.m();
      assert(a.n()==m_);
    }
    if ( transb == 'N' || transb == 'n' )
    {
      n = b.n();
      assert(k==b.m());
    }
    else
    {
      n = b.m();
      assert(k==b.n());
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dgemm(&transa, &transb, &m, &n, &k, &alpha, a.val, &a.lld_,
            b.val, &b.lld_, &beta, val, &lld_);
    }
    else
    {
      int ione=1;
      pdgemm(&transa, &transb, &m, &n, &k, &alpha,
           a.val, &ione, &ione, a.desc_,
           b.val, &ione, &ione, b.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// complex matrix multiplication
// this = alpha*op(A)*op(B)+beta*this
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::gemm(char transa, char transb,
                complex<double> alpha, const ComplexMatrix& a,
                const ComplexMatrix& b, complex<double> beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( ictxt_ == b.ictxt() );

  if ( active() )
  {
    int m, n, k;
    if ( transa == 'N' || transa == 'n' )
    {
      m = a.m();
      k = a.n();
      assert(a.m()==m_);
    }
    else
    {
      m = a.n();
      k = a.m();
      assert(a.n()==m_);
    }
    if ( transb == 'N' || transb == 'n' )
    {
      n = b.n();
      assert(k==b.m());
    }
    else
    {
      n = b.m();
      assert(k==b.n());
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zgemm(&transa, &transb, &m, &n, &k, &alpha, a.val, &a.lld_,
            b.val, &b.lld_, &beta, val, &lld_);
    }
    else
    {
      int ione=1;
      pzgemm(&transa, &transb, &m, &n, &k, &alpha,
           a.val, &ione, &ione, a.desc_,
           b.val, &ione, &ione, b.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// symmetric_matrix * matrix multiplication
// this = beta * this + alpha * a * b
// this = beta * this + alpha * b * a
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::symm(char side, char uplo,
                double alpha, const DoubleMatrix& a,
                const DoubleMatrix& b, double beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( ictxt_ == b.ictxt() );

  if ( active() )
  {
    assert(a.n()==a.m());
    if ( side == 'L' || side == 'l' )
    {
      assert(a.n()==b.m());
    }
    else
    {
      assert(a.m()==b.n());
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dsymm(&side, &uplo, &m_, &n_, &alpha, a.val, &a.lld_,
            b.val, &b.lld_, &beta, val, &lld_);
    }
    else
    {
      int ione=1;
      pdsymm(&side, &uplo, &m_, &n_, &alpha,
           a.val, &ione, &ione, a.desc_,
           b.val, &ione, &ione, b.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// hermitian_matrix * matrix multiplication
// this = beta * this + alpha * a * b
// this = beta * this + alpha * b * a
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::hemm(char side, char uplo,
                complex<double> alpha, const ComplexMatrix& a,
                const ComplexMatrix& b, complex<double> beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( ictxt_ == b.ictxt() );

  if ( active() )
  {
    assert(a.n()==a.m());
    if ( side == 'L' || side == 'l' )
    {
      assert(a.n()==b.m());
    }
    else
    {
      assert(a.m()==b.n());
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zhemm(&side, &uplo, &m_, &n_, &alpha, a.val, &a.lld_,
            b.val, &b.lld_, &beta, val, &lld_);
    }
    else
    {
      int ione=1;
      pzhemm(&side, &uplo, &m_, &n_, &alpha,
           a.val, &ione, &ione, a.desc_,
           b.val, &ione, &ione, b.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// complex_symmetric_matrix * matrix multiplication
// this = beta * this + alpha * a * b
// this = beta * this + alpha * b * a
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::symm(char side, char uplo,
                complex<double> alpha, const ComplexMatrix& a,
                const ComplexMatrix& b, complex<double> beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( ictxt_ == b.ictxt() );

  if ( active() )
  {
    assert(a.n()==a.m());
    if ( side == 'L' || side == 'l' )
    {
      assert(a.n()==b.m());
    }
    else
    {
      assert(a.m()==b.n());
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zsymm(&side, &uplo, &m_, &n_, &alpha, a.val, &a.lld_,
            b.val, &b.lld_, &beta, val, &lld_);
    }
    else
    {
      int ione=1;
      pzsymm(&side, &uplo, &m_, &n_, &alpha,
           a.val, &ione, &ione, a.desc_,
           b.val, &ione, &ione, b.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Compute a matrix-matrix product for a real triangular
// matrix or its transpose.
// *this = alpha op(A) * *this    if side=='l'
// *this = alpha * *this * op(A)  if side=='r'
// where op(A) = A or trans(A)
// alpha is a scalar, *this is an m by n matrix, and A is a unit or non-unit,
// upper- or lower-triangular matrix.
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::trmm(char side, char uplo, char trans, char diag,
                        double alpha, const DoubleMatrix& a)
{
  if ( active() )
  {
    assert(a.m_==a.n_);
    if ( side=='L' || side=='l' )
    {
      assert(a.n_==m_);
    }
    else
    {
      assert(a.n_==n_);
    }
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dtrmm(&side, &uplo, &trans, &diag,
            &m_, &n_, &alpha, a.val, &a.m_, val, &m_);
    }
    else
    {
      int ione=1;
      pdtrmm(&side, &uplo, &trans, &diag, &m_, &n_,
             &alpha, a.val, &ione, &ione, a.desc_,
             val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Solve op(A) * X = alpha * *this  (if side=='l')
// or    X * op(A) = alpha * *this  (if side=='r')
// where op(A) = A or trans(A)
// alpha is a scalar, *this is an m by n matrix, and A is a unit or non-unit,
// upper- or lower-triangular matrix.
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::trsm(char side, char uplo, char trans, char diag,
                        double alpha, const DoubleMatrix& a)
{
  if ( active() )
  {
    assert(a.m_==a.n_);
    if ( side=='L' || side=='l' )
    {
      assert(a.n_==m_);
    }
    else
    {
      assert(a.n_==n_);
    }
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dtrsm(&side, &uplo, &trans, &diag,
            &m_, &n_, &alpha, a.val, &a.m_, val, &m_);
    }
    else
    {
    int ione=1;
    pdtrsm(&side, &uplo, &trans, &diag, &m_, &n_,
           &alpha, a.val, &ione, &ione, a.desc_,
           val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Solve op(A) * X = alpha * *this  (if side=='l')
// or    X * op(A) = alpha * *this  (if side=='r')
// where op(A) = A or trans(A)
// alpha is a scalar, *this is an m by n matrix, and A is a unit or non-unit,
// upper- or lower-triangular matrix.
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::trsm(char side, char uplo, char trans,
  char diag, complex<double> alpha, const ComplexMatrix& a)
{
  if ( active() )
  {
    assert(a.m_==a.n_);
    if ( side=='L' || side=='l' )
    {
      assert(a.n_==m_);
    }
    else
    {
      assert(a.n_==n_);
    }
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      ztrsm(&side, &uplo, &trans, &diag,
            &m_, &n_, &alpha, a.val, &a.m_, val, &m_);
    }
    else
    {
      int ione=1;
      pztrsm(&side, &uplo, &trans, &diag, &m_, &n_,
             &alpha, a.val, &ione, &ione, a.desc_,
             val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Solves a triangular system of the form A * X = B or
// A**T * X = B, where A is a triangular matrix of  order  N,
// and  B  is an N-by-NRHS matrix.
// Output in B.
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::trtrs(char uplo, char trans, char diag,
                         DoubleMatrix& b) const
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dtrtrs(&uplo, &trans, &diag, &m_, &b.n_, val, &m_,
             b.val, &b.m_, &info);
    }
    else
    {
      int ione=1;
      pdtrtrs(&uplo, &trans, &diag, &m_, &b.n_,
      val, &ione, &ione, desc_,
      b.val, &ione, &ione, b.desc_, &info);
    }
    if(info!=0)
    {
      cout <<" Matrix::trtrs, info=" << info << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 2);
#else
      exit(2);
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Solves a triangular system of the form A * X = B or
// A**T * X = B, where A is a triangular matrix of  order  N,
// and  B  is an N-by-NRHS matrix.
// Output in B.
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::trtrs(char uplo, char trans, char diag,
                         ComplexMatrix& b) const
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      ztrtrs(&uplo, &trans, &diag, &m_, &b.n_, val, &m_,
             b.val, &b.m_, &info);
    }
    else
    {
      int ione=1;
      pztrtrs(&uplo, &trans, &diag, &m_, &b.n_,
      val, &ione, &ione, desc_,
      b.val, &ione, &ione, b.desc_, &info);
    }
    if(info!=0)
    {
      cout <<" ComplexMatrix::trtrs, info=" << info << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 2);
#else
      exit(2);
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// LU decomposition of a double matrix
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::lu(valarray<int>& ipiv)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);
    ipiv.resize(mloc_+mb_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dgetrf(&m_, &n_, val, &m_, &ipiv[0], &info);
    }
    else
    {
      int ione=1;
      pdgetrf(&m_, &n_, val, &ione, &ione, desc_, &ipiv[0], &info);
    }
    if(info!=0)
    {
      cout << " DoubleMatrix::lu, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// LU decomposition of a complex matrix
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::lu(valarray<int>& ipiv)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);
    ipiv.resize(mloc_+mb_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zgetrf(&m_, &n_, val, &m_, &ipiv[0], &info);
    }
    else
    {
      int ione=1;
      pzgetrf(&m_, &n_, val, &ione, &ione, desc_, &ipiv[0], &info);
    }
    if(info!=0)
    {
      cout << " ComplexMatrix::lu, info=" << info << endl;
#ifdef USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 2);
#else
      exit(2);
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// inverse of a square double matrix
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::inverse(void)
{
  if ( active() )
  {
    assert(m_==n_);
    valarray<int> ipiv(mloc_+mb_);
    // LU decomposition
    lu(ipiv);
    inverse_from_lu(ipiv);
  }
}

////////////////////////////////////////////////////////////////////////////////
// determinant of a square double matrix in LU form
////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::det_from_lu(valarray<int> ipiv)
{
  if ( active() )
  {
    assert(m_==n_);

    // compute determinant
    valarray<double> diag(n_);
    for ( int ii = 0; ii < n_; ii++ )
    {
      int iii = l(ii) * mb_ + x(ii);
      int jjj = m(ii) * nb_ + y(ii);
      if ( pr(ii) == ctxt_.myrow()
          && pc(ii) == ctxt_.mycol() )
        diag[ii] = val[iii+mloc_*jjj];
    }
    ctxt_.dsum(n_,1,(double*)&diag[0],n_);

    double det = 1.0;
    for ( int ii = 0; ii < n_; ii++ )
      det *= diag[ii];
    det *= signature(ipiv);

    return det;
  }
  return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// inverse and determinant of a square double matrix
////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::inverse_det(void)
{
  if ( active() )
  {
    assert(m_==n_);
    valarray<int> ipiv(mloc_+mb_);
    lu(ipiv);
    double det = det_from_lu(ipiv);
    inverse_from_lu(ipiv);
    return det;
  }
  return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// inverse from an LU decomposed square matrix
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::inverse_from_lu(valarray<int>& ipiv)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    // Compute inverse using LU decomposition and array ipiv
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      valarray<double> work(1);
      int lwork = -1;
      // First call to compute optimal size of work array, returned in work[0]
      dgetri(&m_, val, &m_, &ipiv[0], &work[0], &lwork, &info);
      lwork = (int) work[0] + 1;
      work.resize(lwork);
      dgetri(&m_, val, &m_, &ipiv[0], &work[0], &lwork, &info);
    }
    else
    {
      valarray<double> work(1);
      valarray<int> iwork(1);
      int lwork = -1;
      int liwork = -1;
      int ione = 1;
      // First call to compute dimensions of work arrays lwork and liwork
      // dimensions are returned in work[0] and iwork[0];
      pdgetri(&n_, val, &ione, &ione, desc_, &ipiv[0],
              &work[0], &lwork, &iwork[0], &liwork, &info);
      lwork = (int) work[0] + 1;
      liwork = iwork[0];
      work.resize(lwork);
      iwork.resize(liwork);

      // Compute inverse
      pdgetri(&n_, val, &ione, &ione, desc_, &ipiv[0],
              &work[0], &lwork, &iwork[0], &liwork, &info);
    }
    if(info!=0)
    {
      cout << " DoubleMatrix::inverse_from_lu, info(getri)=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// inverse of a complex square matrix
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::inverse(void)
{
  valarray<int> ipiv;
  lu(ipiv);
  inverse_from_lu(ipiv);
}

////////////////////////////////////////////////////////////////////////////////
// determinant of a complex square matrix in LU form
////////////////////////////////////////////////////////////////////////////////
complex<double> ComplexMatrix::det_from_lu(valarray<int> ipiv)
{
  if ( active() )
  {
    assert(m_==n_);

    // compute determinant
    valarray<complex<double> > diag(n_);
    for ( int ii = 0; ii < n_; ii++ )
    {
      int iii = l(ii) * mb_ + x(ii);
      int jjj = m(ii) * nb_ + y(ii);
      if ( pr(ii) == ctxt_.myrow() && pc(ii) == ctxt_.mycol() )
        diag[ii] = val[iii+mloc_*jjj];
    }
    ctxt_.dsum(n_*2,1,(double*)&diag[0],n_*2);

    complex<double> det = 1.0;
    for ( int ii = 0; ii < n_; ii++ )
      det *= diag[ii];
    det *= signature(ipiv);

    return det;
  }
  return complex<double>(0.0,0.0);
}

////////////////////////////////////////////////////////////////////////////////
// inverse and determinant of a complex square matrix
////////////////////////////////////////////////////////////////////////////////
complex<double> ComplexMatrix::inverse_det(void)
{
  if ( active() )
  {
    assert(m_==n_);
    valarray<int> ipiv(mloc_+mb_);
    lu(ipiv);
    complex<double> det = det_from_lu(ipiv);
    inverse_from_lu(ipiv);
    return det;
  }
  return complex<double>(0.0,0.0);
}

////////////////////////////////////////////////////////////////////////////////
// compute inverse of an LU decomposed matrix
void ComplexMatrix::inverse_from_lu(valarray<int>& ipiv)
{
  // it is assumed that the current matrix is LU decomposed
  int info;
  if ( active() )
  {
    assert(m_==n_);
    // Compute inverse using LU decomposition and array ipiv computed in lu()
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      valarray< complex<double> > work(1);
      int lwork = -1;
      // First call to compute optimal size of work array, returned in work[0]
      zgetri(&m_, val, &m_, &ipiv[0], &work[0], &lwork, &info);
      lwork = (int) work[0].real() + 1;
      work.resize(lwork);
      zgetri(&m_, val, &m_, &ipiv[0], &work[0], &lwork, &info);
    }
    else
    {
      valarray< complex<double> > work(1);
      valarray<int> iwork(1);
      int lwork = -1;
      int liwork = -1;
      int ione = 1;
      pzgetri(&n_, val, &ione, &ione, desc_, &ipiv[0],
              &work[0], &lwork, &iwork[0], &liwork, &info);
      lwork = (int) work[0].real() + 1;
      liwork = iwork[0];
      work.resize(lwork);
      iwork.resize(liwork);

      // Compute inverse
      pzgetri(&n_, val, &ione, &ione, desc_, &ipiv[0],
              &work[0], &lwork, &iwork[0], &liwork, &info);
    }
    if(info!=0)
    {
      cout << " ComplexMatrix::inverse, info(getri)=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Real Cholesky factorization of a
// symmetric positive definite distributed matrix
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::potrf(char uplo)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dpotrf(&uplo, &m_, val, &m_, &info);
    }
    else
    {
      int ione=1;
      pdpotrf(&uplo, &m_, val, &ione, &ione, desc_, &info);
    }
    if(info!=0)
    {
      cout << " DoubleMatrix::potrf, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Complex Cholesky factorization of a
// hermitian positive definite distributed matrix
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::potrf(char uplo)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zpotrf(&uplo, &m_, val, &m_, &info);
    }
    else
    {
      int ione=1;
      pzpotrf(&uplo, &m_, val, &ione, &ione, desc_, &info);
    }
    if(info!=0)
    {
      cout << " ComplexMatrix::potrf, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Compute  the inverse of a real symmetric positive definite matrix
// using the Cholesky factorization A = U**T*U or A = L*L**T computed
// by DoubleMatrix::potrf
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::potri(char uplo)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dpotri(&uplo, &m_, val, &m_, &info);
    }
    else
    {
      int ione=1;
      pdpotri(&uplo, &m_, val, &ione, &ione, desc_, &info);
    }
    if(info!=0)
    {
      cout << " Matrix::potri, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Inverse of a triangular matrix
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::trtri(char uplo, char diag)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dtrtri(&uplo, &diag, &m_, val, &m_, &info);
    }
    else
    {
      int ione=1;
      pdtrtri(&uplo, &diag, &m_, val, &ione, &ione, desc_, &info);
    }
    if(info!=0)
    {
      cout << " Matrix::trtri, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::trtri(char uplo, char diag)
{
  int info;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      ztrtri(&uplo, &diag, &m_, val, &m_, &info);
    }
    else
    {
      int ione=1;
      pztrtri(&uplo, &diag, &m_, val, &ione, &ione, desc_, &info);
    }
    if(info!=0)
    {
      cout << " Matrix::trtri, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Polar decomposition A = UH
// Replace *this with its orthogonal polar factor U
// return when iter > maxiter or ||I - X^T*X|| < tol
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::polar(double tol, int maxiter)
{
  DoubleMatrix x(ctxt_,m_,n_,mb_,nb_);
  DoubleMatrix xp(ctxt_,m_,n_,mb_,nb_);

  DoubleMatrix q(ctxt_,n_,n_,nb_,nb_);
  DoubleMatrix qt(ctxt_,n_,n_,nb_,nb_);
  DoubleMatrix t(ctxt_,n_,n_,nb_,nb_);

  double qnrm2 = numeric_limits<double>::max();
  int iter = 0;
  x = *this;
  while ( iter < maxiter && qnrm2 > tol )
  {
    // q = I - x^T * x
    q.identity();
    q.syrk('l','t',-1.0,x,1.0);
    q.symmetrize('l');

    double qnrm2 = q.nrm2();
#ifdef DEBUG
    if ( ctxt_.onpe0() )
      cout << " DoubleMatrix::polar: qnrm2 = " << qnrm2 << endl;
#endif

    // choose Bjork-Bowie or Higham iteration depending on q.nrm2

    // threshold value
    // see A. Bjork and C. Bowie, SIAM J. Num. Anal. 8, 358 (1971) p.363
    if ( qnrm2 < 1.0 )
    {
      // Bjork-Bowie iteration
      // compute xp = x * ( I + 0.5*q * ( I + 0.75 * q ) )

      // t = ( I + 0.75 * q )
      t.identity();
      t.axpy(0.75,q);

      // compute q*t
      qt.gemm('n','n',1.0,q,t,0.0);

      // xp = x * ( I + 0.5*q * ( I + 0.75 * q ) )
      //    = x * ( I + 0.5 * qt )
      // Use t to store (I + 0.5 * qt)
      t.identity();
      t.axpy(0.5,qt);

      // t now contains (I + 0.5 * qt)
      // xp = x * t
      xp.gemm('n','n',1.0,x,t,0.0);

      // update x
      x = xp;
    }
    else
    {
      // Higham iteration
      assert(m_==n_);
      //if ( ctxt_.onpe0() )
      //  cout << " DoubleMatrix::polar: using Higham algorithm" << endl;
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();
      // t now contains X^-T
      // xp = 0.5 * ( x + x^-T );
      for ( int i = 0; i < x.size(); i++ )
        x[i] = 0.5 * ( x[i] + t[i] );
    }
    iter++;
  }
  *this = x;
}

////////////////////////////////////////////////////////////////////////////////
// Polar decomposition A = UH (complex case)
// Replace *this with its unitary polar factor U
// return when iter > maxiter or ||I - X^H*X|| < tol
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::polar(double tol, int maxiter)
{
  ComplexMatrix x(ctxt_,m_,n_,mb_,nb_);
  ComplexMatrix xp(ctxt_,m_,n_,mb_,nb_);

  ComplexMatrix q(ctxt_,n_,n_,nb_,nb_);
  ComplexMatrix qt(ctxt_,n_,n_,nb_,nb_);
  ComplexMatrix t(ctxt_,n_,n_,nb_,nb_);

  double qnrm2 = numeric_limits<double>::max();
  int iter = 0;
  x = *this;
  while ( iter < maxiter && qnrm2 > tol )
  {
    // q = I - x^T * x
    q.identity();
    q.herk('l','c',-1.0,x,1.0);
    q.symmetrize('l');

    double qnrm2 = q.nrm2();
#ifdef DEBUG
    if ( ctxt_.onpe0() )
      cout << " ComplexMatrix::polar: qnrm2 = " << qnrm2 << endl;
#endif

    // choose Bjork-Bowie or Higham iteration depending on q.nrm2

    // threshold value
    // see A. Bjork and C. Bowie, SIAM J. Num. Anal. 8, 358 (1971) p.363
    if ( qnrm2 < 1.0 )
    {
      // Bjork-Bowie iteration
      // compute xp = x * ( I + 0.5*q * ( I + 0.75 * q ) )

      // t = ( I + 0.75 * q )
      t.identity();
      t.axpy(0.75,q);

      // compute q*t
      qt.gemm('n','n',1.0,q,t,0.0);

      // xp = x * ( I + 0.5*q * ( I + 0.75 * q ) )
      //    = x * ( I + 0.5 * qt )
      // Use t to store (I + 0.5 * qt)
      t.identity();
      t.axpy(0.5,qt);

      // t now contains (I + 0.5 * qt)
      // xp = x * t
      xp.gemm('n','n',1.0,x,t,0.0);

      // update x
      x = xp;
    }
    else
    {
      // Higham iteration
      assert(m_==n_);
      //if ( ctxt_.onpe0() )
      //  cout << " ComplexMatrix::polar: using Higham algorithm" << endl;
      // t = X^H
      t.transpose(1.0,x,0.0);
      t.inverse();
      // t now contains X^-H
      // xp = 0.5 * ( x + x^-H );
      for ( int i = 0; i < x.size(); i++ )
        x[i] = 0.5 * ( x[i] + t[i] );
    }
    iter++;
  }
  *this = x;
}

////////////////////////////////////////////////////////////////////////////////
// estimate the reciprocal of the condition number (in the 1-norm) of a
// real symmetric positive definite matrix using the Cholesky factorization
// A  = U**T*U or A = L*L**T computed by DoubleMatrix::potrf
////////////////////////////////////////////////////////////////////////////////
double DoubleMatrix::pocon(char uplo) const
{
  int info;
  double  rcond=1.;
  double  anorm=1.;
  if ( active() )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      double* work=new double[3*m_];
      int*    iwork=new int[m_];
      dpocon(&uplo, &m_, val, &m_, &anorm, &rcond, work, iwork, &info);
      delete[] iwork;
      delete[] work;
    }
    else
    {
      int     ione=1;
      int     lwork=2*mloc_+3*nloc_+nb_;
      int     liwork=mloc_;
      double* work=new double[lwork];
      int*    iwork=new int[liwork];
      pdpocon(&uplo, &m_, val, &ione, &ione, desc_,
              &anorm, &rcond, work, &lwork, iwork, &liwork, &info);
      if (info!=0)
      {
          cout << "DoubleMatrix::pocon: lwork=" << lwork
               << ", but should be at least " << work[0] << endl;
          cout << "DoubleMatrix::pocon: liwork=" << liwork
               << ", but should be at least " << iwork[0] << endl;
      }
      delete[] iwork;
      delete[] work;
    }
    if(info!=0)
    {
      cout << " Matrix::pocon, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
  return rcond;
}

////////////////////////////////////////////////////////////////////////////////
// symmetric rank k update
// this = beta * this + alpha * A * A^T  (trans=='n')
// this = beta * this + alpha * A^T * A  (trans=='t')
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::syrk(char uplo, char trans,
                double alpha, const DoubleMatrix& a, double beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( n_ == m_ ); // *this must be a square matrix

  if ( active() )
  {
    int n, k;
    if ( trans == 'N' || trans == 'n' )
    {
      n = m_;
      k = a.n();
    }
    else
    {
      n = m_;
      k = a.m();
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dsyrk(&uplo, &trans, &n, &k, &alpha, a.val, &a.m_, &beta, val, &m_);
    }
    else
    {
      int ione = 1;
      pdsyrk(&uplo, &trans, &n, &k, &alpha,
           a.val, &ione, &ione, a.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// hermitian rank k update
// this = beta * this + alpha * A * A^H  (trans=='n')
// this = beta * this + alpha * A^H * A  (trans=='c')
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::herk(char uplo, char trans,
  double alpha, const ComplexMatrix& a, double beta)
{
  assert( ictxt_ == a.ictxt() );
  assert( n_ == m_ ); // *this must be a square matrix

  if ( active() )
  {
    int n, k;
    if ( trans == 'N' || trans == 'n' )
    {
      n = m_;
      k = a.n();
    }
    else if ( trans == 'C' || trans == 'c' )
    {
      n = m_;
      k = a.m();
    }
    else
    {
      cout << " Matrix::herk: invalid parameter trans" << endl;
      MPI_Abort(MPI_COMM_WORLD,2);
    }

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      zherk(&uplo, &trans, &n, &k, &alpha, a.val, &a.m_,
            &beta, val, &m_);
    }
    else
    {
      int ione = 1;
      pzherk(&uplo, &trans, &n, &k, &alpha,
           a.val, &ione, &ione, a.desc_,
           &beta, val, &ione, &ione, desc_);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Generate a duplicated matrix from a distributed matrix
//
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::matgather(double *a, int lda) const
{
  if ( active_ )
  {
    memset(a,0,lda*n_*sizeof(double));

    if ( active_ )
    {
      for ( int li=0; li<mblocks(); li++)
      {
        for ( int lj=0; lj<nblocks(); lj++)
        {
          for ( int ii=0; ii<mbs(li);  ii++)
          {
            for ( int jj=0; jj<nbs(lj); jj++)
            {
              assert(i(li,ii)<lda);
              assert((ii+li*mb_)<mloc_);
              assert((jj+lj*nb_)<nloc_);
              a[ i(li,ii) + j(lj,jj)*lda ]
                  = val[ (ii+li*mb_)+(jj+lj*nb_)*mloc_ ];
            }
          }
        }
      }
    }

    int     max_size=200000;
    int     size    = lda*n_;
    int     ione=1;
    int     nblocks = size / max_size;
    int     sizer   = (size % max_size);
    double* work = new double[max_size];

    double  *ptra = a;
    for ( int steps = 0; steps < nblocks; steps++ )
    {
      MPI_Allreduce(ptra, work, max_size, MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
      dcopy(&max_size, work, &ione, ptra, &ione);
      ptra += max_size;
    }
    if ( sizer != 0 )
    {
      MPI_Allreduce(ptra, work, sizer, MPI_DOUBLE, MPI_SUM, ctxt_.comm() );
      dcopy(&sizer, work, &ione, ptra, &ione);
    }

    delete[] work;
  }
}


////////////////////////////////////////////////////////////////////////////////
// initdiag: initialize diagonal elements using a replicated array a[i]
void DoubleMatrix::initdiag(const double* const dmat)
{
  if ( active() )
  {
    // initialize diagonal elements
    if ( active() )
    {
      // loop through all local blocks (ll,mm)
      for ( int ll = 0; ll < mblocks(); ll++)
      {
        for ( int mm = 0; mm < nblocks(); mm++)
        {
          // check if block (ll,mm) has diagonal elements
          int imin = i(ll,0);
          int imax = imin + mbs(ll)-1;
          int jmin = j(mm,0);
          int jmax = jmin + nbs(mm)-1;
          // cout << " process (" << myrow_ << "," << mycol_ << ")"
          // << " block (" << ll << "," << mm << ")"
          // << " imin/imax=" << imin << "/" << imax
          // << " jmin/jmax=" << jmin << "/" << jmax << endl;

          if ((imin <= jmax) && (imax >= jmin))
          {
            // block (ll,mm) holds diagonal elements
            int idiagmin = max(imin,jmin);
            int idiagmax = min(imax,jmax);

            // cout << " process (" << myrow_ << "," << mycol_ << ")"
            // << " holds diagonal elements " << idiagmin << " to " <<
            // idiagmax << " in block (" << ll << "," << mm << ")" << endl;

            for ( int ii = idiagmin; ii <= idiagmax; ii++ )
            {
              // access element (ii,ii)
              int jj = ii;
              int iii = ll * mb_ + x(ii);
              int jjj = mm * nb_ + y(jj);
              val[iii+mloc_*jjj] = dmat[ii];
            }
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// trace
double DoubleMatrix::trace(void) const
{
  assert(m_==n_);

  if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
  {
    double trace = 0.0;
    for ( int i = 0; i < n_; i++ )
      trace += val[i*m_];
    return trace;
  }
  else
  {
    int ione=1;
    return pdlatra(&n_,val,&ione,&ione,desc_);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Reduces a real  symmetric-definite  generalized  eigenproblem  to  standard
// form.  If itype = 1, the problem is A*x = lambda*B*x,
// and A (=*this) is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
// If itype = 2 or 3, the problem is A*B*x = lambda*x or
// B*A*x = lambda*x, and *this is overwritten by U*A*U**T or L**T*A*L.
// B must have been previously factorized as U**T*U or L*L**T by
// DoubleMatrix::dpotrf.
////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::sygst(int itype, char uplo, const DoubleMatrix& b)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      dsygst(&itype, &uplo, &m_, val, &m_, b.val, &b.m_, &info);
    }
    else
    {
      int ione=1;
      double  scale;
      pdsygst(&itype, &uplo, &m_, val, &ione, &ione, desc_,
      b.val, &ione, &ione, b.desc_, &scale, &info);
    }
    if ( info != 0 )
    {
      cout << " Matrix::sygst, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues and eigenvectors of *this
// store eigenvalues in w, eigenvectors in z
void DoubleMatrix::syev(char uplo, valarray<double>& w, DoubleMatrix& z)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'V';
    double* work;
    int lwork;
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      lwork=-1;
      double tmplwork;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work = new double[lwork];
      z = *this;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], work, &lwork, &info);
    }
    else
    {
      int ione=1;
      lwork=-1;
      double tmpwork;
      pdsyev(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, &tmpwork, &lwork,
              &info);
      lwork = (int) (tmpwork + 1);
      // set lwork to max value among all tasks
      ctxt_.imax(1,1,&lwork,1);
      work=new double[lwork];
      pdsyev(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, work, &lwork,
              &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
    }
    if ( info != 0 )
    {
      cout << " Matrix::syev requires lwork>=" << work[0] << endl;
      cout << " Matrix::syev, lwork>=" << lwork << endl;
      cout << " Matrix::syev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues and eigenvectors of *this
// store eigenvalues in w, eigenvectors in z
// using the divide and conquer algorithm of Tisseur and Dongarra
void DoubleMatrix::syevd(char uplo, valarray<double>& w, DoubleMatrix& z)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'V';
    int lwork;
    double* work;
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      lwork=-1;
      double tmplwork;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work = new double[lwork];
      z = *this;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], work, &lwork, &info);
    }
    else
    {
      int ione=1;
      lwork=-1;
      double tmpwork;
      int liwork=-1;
      int tmpiwork;
      pdsyevd(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, &tmpwork, &lwork,
              &tmpiwork, &liwork, &info);
      lwork = (int) (tmpwork + 1);
      // set lwork to max value among all tasks
      ctxt_.imax(1,1,&lwork,1);
      work=new double[lwork];
      liwork = tmpiwork;
      int* iwork = new int[liwork];
      pdsyevd(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
            z.val, &ione, &ione, z.desc_, work, &lwork, iwork, &liwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      delete[] iwork;
    }
    if ( info != 0 )
    {
      cout << " Matrix::syev requires lwork>=" << work[0] << endl;
      cout << " Matrix::syev, lwork>=" << lwork << endl;
      cout << " Matrix::syev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues and eigenvectors of *this
// store eigenvalues in w, eigenvectors in z
// using the expert driver
void DoubleMatrix::syevx(char uplo, valarray<double>& w, DoubleMatrix& z,
                         double abstol)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'V';
    int lwork;
    double* work;
    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      lwork=-1;
      double tmplwork;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work = new double[lwork];
      z = *this;
      dsyev(&jobz, &uplo, &m_, z.val, &m_, &w[0], work, &lwork, &info);
    }
    else
    {
      char range = 'A';
      int ione=1;
      lwork=-1;
      double tmpwork;
      int liwork=-1;
      int tmpiwork;
      valarray<int> ifail(n_);
      int nfound=-1;
      int nz=-1;
      int il=1, iu=n_;
      double vl=0, vu=0;
      double orfac=-1.0;
      valarray<int> icluster(2*ctxt_.size());
      valarray<double> gap(ctxt_.size());
      pdsyevx(&jobz, &range, &uplo, &m_, val, &ione, &ione, desc_,
              &vl, &vu, &il, &iu, &abstol, &nfound, &nz, &w[0],
              &orfac, z.val, &ione, &ione, z.desc_, &tmpwork, &lwork,
              &tmpiwork, &liwork, &ifail[0], &icluster[0], &gap[0], &info);
      assert(info==0);
      lwork = (int) (tmpwork + 1);
      work=new double[lwork];
      liwork = tmpiwork;
      int* iwork = new int[liwork];
      pdsyevx(&jobz, &range, &uplo, &m_, val, &ione, &ione, desc_,
              &vl, &vu, &il, &iu, &abstol, &nfound, &nz, &w[0],
              &orfac, z.val, &ione, &ione, z.desc_, work, &lwork,
              iwork, &liwork, &ifail[0], &icluster[0], &gap[0], &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      delete[] iwork;
    }
    if ( info != 0 )
    {
      cout << " Matrix::syev requires lwork>=" << work[0] << endl;
      cout << " Matrix::syev, lwork>=" << lwork << endl;
      cout << " Matrix::syev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues (only) of *this
// store eigenvalues in w
void DoubleMatrix::syev(char uplo, valarray<double>& w)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'N';
    int lwork;
    double* work;

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      lwork=-1;
      double tmplwork;
      dsyev(&jobz, &uplo, &m_, val, &m_, &w[0], &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work = new double[lwork];
      dsyev(&jobz, &uplo, &m_, val, &m_, &w[0], work, &lwork, &info);
    }
    else
    {
      int ione=1;
      lwork=-1;
      double tmplwork;
      double *zval = 0; // zval is not referenced since jobz == 'N'
      int * descz = 0;
      pdsyev(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work=new double[lwork];
      pdsyev(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, work, &lwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
    }
    if ( info != 0 )
    {
      cout << " Matrix::syev requires lwork>=" << work[0] << endl;
      cout << " Matrix::syev, lwork>=" << lwork << endl;
      cout << " Matrix::syev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues (only) of *this
// store eigenvalues in w
// using the divide and conquer method of Tisseur and Dongarra
void DoubleMatrix::syevd(char uplo, valarray<double>& w)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'N';
    int lwork,liwork;
    double* work;

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      lwork=-1;
      double tmplwork;
      dsyev(&jobz, &uplo, &m_, val, &m_, &w[0], &tmplwork, &lwork, &info);
      lwork = (int) tmplwork + 1;
      work = new double[lwork];
      dsyev(&jobz, &uplo, &m_, val, &m_, &w[0], work, &lwork, &info);
    }
    else
    {
      int ione=1;
      lwork=-1;
      double tmpwork;
      liwork=-1;
      int tmpiwork;
      double *zval = 0; // zval is not referenced since jobz == 'N'
      int * descz = 0;
      pdsyevd(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
              zval, &ione, &ione, descz, &tmpwork, &lwork,
              &tmpiwork, &liwork, &info);
      lwork = (int) (tmpwork + 1);
      work=new double[lwork];
      liwork = tmpiwork;
      int* iwork = new int[liwork];
      pdsyevd(&jobz, &uplo, &m_, val, &ione, &ione, desc_, &w[0],
              zval, &ione, &ione, descz, work, &lwork, iwork, &liwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      delete[] iwork;
    }
    if ( info != 0 )
    {
      cout << " Matrix::syev requires lwork>=" << work[0] << endl;
      cout << " Matrix::syev, lwork>=" << lwork << endl;
      cout << " Matrix::syev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::heev(char uplo, valarray<double>& w, ComplexMatrix& z)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'V';
    int lwork;
    complex<double>* work;
    double* rwork;

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      // request optimal lwork size
      int lwork=-1;
      complex<double> tmplwork;
      int lrwork = max(1,3*n_-2);
      rwork = new double[lrwork];
      zheev(&jobz, &uplo, &m_, z.val, &m_, &w[0], &tmplwork, &lwork,
            rwork, &info);
      lwork = (int) real(tmplwork) + 1;
      work = new complex<double>[lwork];
      z=*this;
      zheev(&jobz, &uplo, &m_, z.val, &m_, &w[0], work, &lwork,
            rwork, &info);
    }
    else
    {
      int ione=1;
      lwork=-1;
      int lrwork=-1;
      complex<double> tmplwork;
      double tmplrwork;
      // first call to get optimal lwork and lrwork sizes
      pzheev(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, &tmplwork, &lwork,
              &tmplrwork, &lrwork, &info);
      lwork = (int) real(tmplwork) + 1;
      work = new complex<double>[lwork];
      lrwork = (int) tmplrwork + 1;
      // direct calculation of lrwork to avoid bug in pzheev
      lrwork = 1 + 9*n_ + 3*mloc_*nloc_;
      rwork = new double[lrwork];
      pzheev(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, work, &lwork,
              rwork, &lrwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
    }
    if ( info != 0 )
    {
      cout << " Matrix::heev requires lwork>=" << work[0] << endl;
      cout << " Matrix::heev, lwork>=" << lwork << endl;
      cout << " Matrix::heev, info=" << info<< endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
    delete[] rwork;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::heevd(char uplo, valarray<double>& w, ComplexMatrix& z)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'V';
    int lwork,liwork;
    complex<double>* work;
    double* rwork;

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      // request optimal lwork size
      lwork=-1;
      complex<double> tmplwork;
      int lrwork = max(1,3*n_-2);
      rwork = new double[lrwork];
      zheev(&jobz, &uplo, &m_, z.val, &m_, &w[0], &tmplwork, &lwork,
            rwork, &info);
      lwork = (int) real(tmplwork) + 1;
      work = new complex<double>[lwork];
      z=*this;
      zheev(&jobz, &uplo, &m_, z.val, &m_, &w[0], work, &lwork,
            rwork, &info);
      if ( info != 0 )
      {
        cout << " Matrix::heevd requires lwork>=" << work[0] << endl;
        cout << " Matrix::heevd, lwork>=" << lwork << endl;
        cout << " Matrix::heevd, info=" << info<< endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
      }
    }
    else
    {
      int ione=1;
      lwork=-1;
      int lrwork=-1;
      liwork=-1;
      complex<double> tmplwork;
      double tmplrwork;
      int tmpliwork;
      // first call to get optimal lwork and lrwork sizes
      pzheevd(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, &tmplwork, &lwork,
              &tmplrwork, &lrwork, &tmpliwork, &liwork, &info);
      lwork = (int) real(tmplwork) + 1;
      work = new complex<double>[lwork];
      lrwork = (int) tmplrwork + 1;
      rwork = new double[lrwork];
      liwork = tmpliwork;
      int* iwork = new int[liwork];
      pzheevd(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
              z.val, &ione, &ione, z.desc_, work, &lwork,
              rwork, &lrwork, iwork, &liwork, &info);
      //MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      if ( info != 0 )
      {
        cout << " Matrix::heevd requires lwork>=" << work[0] << endl;
        cout << " Matrix::heevd, lwork>=" << lwork << endl;
        cout << " Matrix::heevd, liwork>=" << liwork << endl;
        cout << " Matrix::heevd, info=" << info<< endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
        delete[] work;
        delete[] rwork;
        delete[] iwork;
      }
    }
    delete[] work;
    delete[] rwork;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues (only) of hermitian matrix *this
void ComplexMatrix::heev(char uplo, valarray<double>& w)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'N';
    int lwork;
    complex<double>* work;

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      // request optimal lwork size
      lwork=-1;
      complex<double> tmplwork;
      int lrwork = max(1,3*n_-2);
      double* rwork = new double[lrwork];
      zheev(&jobz, &uplo, &m_, val, &m_, &w[0], &tmplwork, &lwork,
            rwork, &info);
      lwork = (int) real(tmplwork);
      work = new complex<double>[lwork];
      zheev(&jobz, &uplo, &m_, val, &m_, &w[0], work, &lwork,
            rwork, &info);
      delete[] rwork;
    }
    else
    {
      int ione=1;
      lwork=-1;
      int lrwork=-1;
      complex<double> tmplwork;
      double tmplrwork;
      complex<double> *zval = 0;
      int *descz = 0;
      // first call to get optimal lwork and lrwork sizes
      pzheev(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, &tmplwork, &lwork,
             &tmplrwork, &lrwork, &info);
      lwork = (int) real(tmplwork) + 1;
      work = new complex<double>[lwork];
      lrwork = (int) tmplrwork + 1;
      double* rwork = new double[lrwork];
      pzheev(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, work, &lwork,
             rwork, &lrwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      delete[] rwork;
    }
    if ( info != 0 )
    {
      cout << " Matrix::heev requires lwork>=" << work[0] << endl;
      cout << " Matrix::heev, lwork>=" << lwork << endl;
      cout << " Matrix::heev, info=" << info << endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
    delete[] work;
  }
}

////////////////////////////////////////////////////////////////////////////////
// compute eigenvalues (only) of hermitian matrix *this
void ComplexMatrix::heevd(char uplo, valarray<double>& w)
{
  int info;
  if ( active_ )
  {
    assert(m_==n_);
    char jobz = 'N';

    if ( ( nprow_ == 1 ) && ( npcol_ == 1 ) )
    {
      // request optimal lwork size
      int lwork=-1;
      complex<double> tmplwork;
      int lrwork = max(1,3*n_-2);
      double* rwork = new double[lrwork];
      zheev(&jobz, &uplo, &m_, val, &m_, &w[0], &tmplwork, &lwork,
            rwork, &info);
      lwork = (int) real(tmplwork);
      complex<double>* work = new complex<double>[lwork];
      zheev(&jobz, &uplo, &m_, val, &m_, &w[0], work, &lwork,
            rwork, &info);
      if ( info != 0 )
      {
        cout << " Matrix::heevd requires lwork>=" << work[0] << endl;
        cout << " Matrix::heevd, lwork>=" << lwork << endl;
        cout << " Matrix::heevd, lrwork>=" << lrwork << endl;
        cout << " Matrix::heevd, info=" << info << endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
      }
      delete[] work;
      delete[] rwork;
    }
    else
    {
      int ione=1;
      int lwork=-1;
      int lrwork=-1;
      int liwork=-1;
      complex<double> tmplwork;
      double tmplrwork;
      int tmpliwork;
      complex<double> *zval = 0;
      int *descz = 0;
      // first call to get optimal lwork and lrwork sizes
      pzheevd(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, &tmplwork, &lwork,
             &tmplrwork, &lrwork, &tmpliwork, &liwork, &info);
      lwork = (int) real(tmplwork) + 1;
      complex<double>* work = new complex<double>[lwork];
      lrwork = (int) tmplrwork + 1;
      liwork = tmpliwork;
      double* rwork = new double[lrwork];
      int* iwork = new int[liwork];
      pzheevd(&jobz, &uplo, &n_, val, &ione, &ione, desc_, &w[0],
             zval, &ione, &ione, descz, work, &lwork,
             rwork, &lrwork, iwork, &liwork, &info);
      MPI_Bcast(&w[0], m_, MPI_DOUBLE, 0, ctxt_.comm());
      if ( info != 0 )
      {
        cout << " Matrix::heevd requires lwork>=" << work[0] << endl;
        cout << " Matrix::heevd, lwork>=" << lwork << endl;
        cout << " Matrix::heevd, lrwork>=" << lrwork << endl;
        cout << " Matrix::heevd, liwork>=" << liwork << endl;
        cout << " Matrix::heevd, info=" << info << endl;
        MPI_Abort(MPI_COMM_WORLD, 2);
      }
      delete[] work;
      delete[] rwork;
      delete[] iwork;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::lapiv(char direc, char rowcol, const int *ipiv)
{
  // Perform a permutation of rows or columns of *this
  //
  // Permutation of rows: (rowcol='R' or 'r')
  // the array ipiv is distributed over a process column
  // and is replicated on all process columns
  // ipiv has size mloc and contains the local values of the permutation
  //
  // Permutation of columns: (rowcol='C' or 'c')
  // the array ipiv is distributed over a process row
  // and is replicated on all process rows
  // ipiv has size nloc and contains the local values of the permutation

  const bool rowcol_r = ( rowcol=='R' || rowcol=='r' );
  const bool rowcol_c = ( rowcol=='C' || rowcol=='c' );
  assert(rowcol_r || rowcol_c);

  // ipivtmp: extended permutation array for use in pdlapv2
  // (see scalapack documentation)
  vector<int> ipivtmp;
  // descriptor of the ipivtmp distributed vector
  int desc_ip[9];
  if ( rowcol_r )
  {
    // permuting rows
    ipivtmp.resize(mloc_ + mb_);
    for ( int i = 0; i < mloc_; i++)
      ipivtmp[i] = ipiv[i];
    // initialize descriptor: ipivtmp is (mx1)
    desc_ip[0] = 1;     // dtype
    desc_ip[1] = ictxt_;// ctxt
    desc_ip[2] = m_+mb_*nprow_; // m (see details in pldapv2.f)
    desc_ip[3] = 1;     // n
    desc_ip[4] = mb_;   // mb
    desc_ip[5] = 1;     // nb
    desc_ip[6] = 0;     // rsrc
    desc_ip[7] = 0;     // csrc
    desc_ip[8] = mloc_; // lld
  }
  else
  {
    // permuting columns
    ipivtmp.resize(nloc_ + nb_);
    for ( int i = 0; i < nloc_; i++)
      ipivtmp[i] = ipiv[i];
    // initialize descriptor: ipivtmp is (1xn)
    desc_ip[0] = 1;     // dtype
    desc_ip[1] = ictxt_;// ctxt
    desc_ip[2] = 1;     // m
    desc_ip[3] = n_+nb_*npcol_; // n (see details in pdlapv2.f)
    desc_ip[4] = 1;     // mb
    desc_ip[5] = nb_;   // nb
    desc_ip[6] = 0;     // rsrc
    desc_ip[7] = 0;     // csrc
    desc_ip[8] = 1;     // lld
  }

  int one = 1;
  pdlapv2(&direc, &rowcol, &m_, &n_, val, &one, &one, desc_,
          &ipivtmp[0], &one, &one, desc_ip);
}

////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::lapiv(char direc, char rowcol, const int *ipiv)
{
  // Perform a permutation of rows or columns of *this
  //
  // Permutation of rows: (rowcol='R' or 'r')
  // the array ipiv is distributed over a process column
  // and is replicated on all process columns
  // ipiv has size mloc and contains the local values of the permutation
  //
  // Permutation of columns: (rowcol='C' or 'c')
  // the array ipiv is distributed over a process row
  // and is replicated on all process rows
  // ipiv has size nloc and contains the local values of the permutation

  const bool rowcol_r = ( rowcol=='R' || rowcol=='r' );
  const bool rowcol_c = ( rowcol=='C' || rowcol=='c' );
  assert(rowcol_r || rowcol_c);

  // ipivtmp: extended permutation array for use in pdlapv2
  // (see scalapack documentation)
  vector<int> ipivtmp;
  // descriptor of the ipivtmp distributed vector
  int desc_ip[9];
  if ( rowcol_r )
  {
    // permuting rows
    ipivtmp.resize(mloc_ + mb_);
    for ( int i = 0; i < mloc_; i++)
      ipivtmp[i] = ipiv[i];
    // initialize descriptor: ipivtmp is (mx1)
    desc_ip[0] = 1;     // dtype
    desc_ip[1] = ictxt_;// ctxt
    desc_ip[2] = m_+mb_*nprow_; // m (see details in pldapv2.f)
    desc_ip[3] = 1;     // n
    desc_ip[4] = mb_;   // mb
    desc_ip[5] = 1;     // nb
    desc_ip[6] = 0;     // rsrc
    desc_ip[7] = 0;     // csrc
    desc_ip[8] = mloc_; // lld
  }
  else
  {
    // permuting columns
    ipivtmp.resize(nloc_ + nb_);
    for ( int i = 0; i < nloc_; i++)
      ipivtmp[i] = ipiv[i];
    // initialize descriptor: ipivtmp is (1x(n+nb))
    desc_ip[0] = 1;     // dtype
    desc_ip[1] = ictxt_;// ctxt
    desc_ip[2] = 1;     // m
    desc_ip[3] = n_+nb_*npcol_; // n (see details in pdlapv2.f)
    desc_ip[4] = 1;     // mb
    desc_ip[5] = nb_;   // nb
    desc_ip[6] = 0;     // rsrc
    desc_ip[7] = 0;     // csrc
    desc_ip[8] = 1;     // lld
  }

  int one = 1;
  pzlapv2(&direc, &rowcol, &m_, &n_, val, &one, &one, desc_,
          &ipivtmp[0], &one, &one, desc_ip);
}

////////////////////////////////////////////////////////////////////////////////
void DoubleMatrix::print(ostream& os) const
{
  // Copy blocks of <blocksize> columns and print them on process (0,0)
  if ( m_ == 0 || n_ == 0 ) return;
  Context ctxtl(MPI_COMM_WORLD,1,1);
  const int blockmemsize = 32768; // maximum memory size of a block in bytes
  // compute maximum block size: must be at least 1
  int maxbs = max(1, (int) ((blockmemsize/sizeof(double))/m_));
  DoubleMatrix t(ctxtl,m_,maxbs);
  int nblocks = n_ / maxbs + ( (n_%maxbs == 0) ? 0 : 1 );
  int ia = 0;
  int ja = 0;
  for ( int jb = 0; jb < nblocks; jb++ )
  {
    int blocksize = ( (jb+1) * maxbs > n_ ) ? n_ % maxbs : maxbs;
    t.getsub(*this,t.m(),blocksize,ia,ja);
    ja += blocksize;
    if ( t.active() )
    {
      // this is done only on pe 0
      for ( int jj = 0; jj < blocksize; jj++ )
      {
        for ( int ii = 0; ii < m_; ii++ )
        {
          os << "(" << ii << "," << jj+jb*maxbs << ")="
             << t.val[ii+t.mloc()*jj] << endl;
        }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ComplexMatrix::print(ostream& os) const
{
  // Copy blocks of <blocksize> columns and print them on process (0,0)
  if ( m_ == 0 || n_ == 0 ) return;
  Context ctxtl(MPI_COMM_WORLD,1,1);
  const int blockmemsize = 32768; // maximum memory size of a block in bytes
  // compute maximum block size: must be at least 1
  int maxbs = max(1, (int) ((blockmemsize/sizeof(complex<double>))/m_));
  ComplexMatrix t(ctxtl,m_,maxbs);
  int nblocks = n_ / maxbs + ( (n_%maxbs == 0) ? 0 : 1 );
  int ia = 0;
  int ja = 0;
  for ( int jb = 0; jb < nblocks; jb++ )
  {
    int blocksize = ( (jb+1) * maxbs > n_ ) ? n_ % maxbs : maxbs;
    t.getsub(*this,t.m(),blocksize,ia,ja);
    ja += blocksize;
    if ( t.active() )
    {
      // this is done only on pe 0
      for ( int jj = 0; jj < blocksize; jj++ )
      {
        for ( int ii = 0; ii < m_; ii++ )
        {
          os << "(" << ii << "," << jj+jb*maxbs << ")="
             << t.val[ii+t.mloc()*jj] << endl;
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const DoubleMatrix& a)
{
  a.print(os);
  return os;
}
////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const ComplexMatrix& a)
{
  a.print(os);
  return os;
}

////////////////////////////////////////////////////////////////////////////////
//
// signature: compute the signature of a row permutation
// defined by a distributed pivot vector ipiv
//
// the vector ipiv is computed by the lu decomposition
//
int DoubleMatrix::signature(valarray<int> ipiv)
{
  // count the number of non-trivial transpositions in the local ipiv vector
  int ntrans = 0;
  for ( int i = 0; i < mloc_; i++ )
  {
    if ( ipiv[i] != iglobal(i) )
      ntrans++;
  }
  // accumulate total number of transpositions
  ctxt_.isum('c',1,1,&ntrans,1);
  return 1 - 2 * ((m_ - ntrans)%2);
}

////////////////////////////////////////////////////////////////////////////////
//
// signature: compute the signature of a row permutation
// defined by a distributed pivot vector ipiv
//
// the vector ipiv is computed by the lu decomposition
//
int ComplexMatrix::signature(valarray<int> ipiv)
{
  // count the number of non-trivial transpositions in the local ipiv vector
  int ntrans = 0;
  for ( int i = 0; i < mloc_; i++ )
  {
    if ( ipiv[i] != iglobal(i) )
      ntrans++;
  }
  // accumulate total number of transpositions
  ctxt_.isum('c',1,1,&ntrans,1);
  return 1 - 2 * ((m_ - ntrans)%2);
}
