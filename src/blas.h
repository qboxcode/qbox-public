////////////////////////////////////////////////////////////////////////////////
//
//  BLAS Header file
//
////////////////////////////////////////////////////////////////////////////////
// $Id: blas.h,v 1.2 2003-11-27 01:16:34 fgygi Exp $

#ifndef BLAS_H
#define BLAS_H

#include <complex>
using namespace std;

// default value for most compilers
#define FTN_LINK extern "C"

#ifdef LINUX
#define FTN_LINK extern "C"
#endif

#ifdef SUNOS
#define FTN_LINK extern "C"
#endif

#ifdef OSF1
#define FTN_LINK extern "C"
#endif

#ifdef HPUX
#define FTN_LINK extern "C"
#define dcopy_ dcopy
#define zcopy_ zcopy
#define daxpy_ daxpy
#define ddot_ ddot
#define dasum_ dasum
#define dsbmv_ dsbmv
#define dgemm_ dgemm
#define dgemv_ dgemv
#define dgesv_ dgesv
#define dscal_ dscal
#define dsyev_ dsyev
#define zdscal_ zdscal
#define idamax_ idamax
#endif

#ifdef AIX
//#define FTN_LINK extern "FORTRAN"
#define FTN_LINK extern "C"
#define dcopy_ dcopy
#define zcopy_ zcopy
#define daxpy_ daxpy
#define ddot_ ddot
#define drot_ drot
#define dasum_ dasum
#define dsbmv_ dsbmv
#define dgemm_ dgemm
#define dgesv_ dgesv
#define dgemv_ dgemv
#define dscal_ dscal
#define dsyev_ dsyev
#define zdscal_ zdscal
#define idamax_ idamax
#define dvea_ dvea
#define dyax_ dyax
#define dnaxpy_ dnaxpy
#define dger_ dger
#endif

#ifndef FTN_LINK
#error "blas.h: undefined platform"
#endif

#ifdef __cplusplus
FTN_LINK {
#endif

void dcopy_(int *n, double *x, int *incx, 
double *y, int *incy );
void zcopy_(int *n, complex<double> *x, int *incx, 
complex<double> *y, int *incy );
void daxpy_(int *n, double *alpha, double *x, int *incx,
double *y, int *incy );
double ddot_(const int *n, const double *a, const int *inca, 
const double *b, const int *incb);
void drot_(int*, double*, int*, double*, int*, double*, double*);
void dgemm_(char *ta, char *tb, int *m, int *n, int *k,
double *alpha, double *a, int *lda, double *b, int *ldb,
double *beta, double *c, int *ldc);
void dgemv_( char *ta, int *m, int *n,
                   double *alpha,  double *a, int *tda,
                   double *x,    int *incx,
                   double *beta,   double *y, int *incy );
 
void dger_(int *,int *, double *, double *, int *,
          double *, int *, double *, int *);
 
void dscal_(int *len, double *alpha, double *x, int *incx);
double dasum_(int *len, double *x, int *incx);
int idamax_(int *len, double *x, int *incx);
void dsyev_(char *c1,char *c2,int *n, 
double *a,int *lda, double *wr,
double *wrk,int *lwrk, int *ierr);
void zdscal_(int *n,double *alpha,complex<double> *x,int *incx);
void dgbmv_(char *trans, int *m, int *n,
int *kl, int *ku, double *alpha, double *a, 
int *lda, double *x, int *incx, double *beta,
double *y, int *incy);
void dsbmv_(char *uplo, int *n, int *k,
double *alpha, double *a, int *lda, double *x, int *incx,
double *beta, double *y, int *incy);
void sspev_(char *vec,char *uplo,int *size,double *ap,
double *wr,double *z,int *n,double *wrk,int *ierr);
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, 
double *b, int *ldb, int *info);
 
void dvea_(int*,double*,int*,double*,int*,double*,int*);
void dyax_(int*,double*,double*,int*,double*,int*);
void dnaxpy_(int*,int*,double*,int*,double*,int*,int*,double*,int*,int*);
 
#ifdef __cplusplus
}
#endif
 
#endif
