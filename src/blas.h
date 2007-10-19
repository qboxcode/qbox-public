////////////////////////////////////////////////////////////////////////////////
//
//  BLAS Header file
//
////////////////////////////////////////////////////////////////////////////////
// $Id: blas.h,v 1.6 2007-10-19 16:24:05 fgygi Exp $

#ifndef BLAS_H
#define BLAS_H

#include <complex>

// default value for most compilers
#define FTN_LINK extern "C"

#ifdef ADD_
#define dcopy  dcopy_
#define zcopy  zcopy_
#define daxpy  daxpy_
#define ddot   ddot_
#define dnrm2  dnrm2_
#define drot   drot_
#define dasum  dasum_
#define dsbmv  dsbmv_
#define dgemm  dgemm_
#define dgesv  dgesv_
#define dgemv  dgemv_
#define dscal  dscal_
#define dsyev  dsyev_
#define zdscal zdscal_
#define idamax idamax_
#define dvea   dvea_
#define dyax   dyax_
#define dnaxpy dnaxpy_
#define dger   dger_
#define zgemm  zgemm_
#endif

#ifdef __cplusplus
FTN_LINK {
#endif

void dcopy(int *n, double *x, int *incx,
double *y, int *incy );
void zcopy(int *n, std::complex<double> *x, int *incx,
std::complex<double> *y, int *incy );
void daxpy(int *n, double *alpha, double *x, int *incx,
double *y, int *incy );
double ddot(const int *n, const double *a, const int *inca,
const double *b, const int *incb);
double dnrm2(const int *n, const double *a, const int *inca);
void drot(int*, double*, int*, double*, int*, double*, double*);
void dgemm(char *ta, char *tb, int *m, int *n, int *k,
  double *alpha, double *a, int *lda, double *b, int *ldb,
  double *beta, double *c, int *ldc);
void zgemm(char *ta, char *tb, int *m, int *n, int *k,
  std::complex<double> *alpha, std::complex<double> *a, int *lda,
  std::complex<double> *b, int *ldb,
  std::complex<double> *beta, std::complex<double> *c, int *ldc);
void dgemv( char *ta, int *m, int *n,
                   double *alpha,  double *a, int *tda,
                   double *x,    int *incx,
                   double *beta,   double *y, int *incy );

void dger(int *,int *, double *, double *, int *,
          double *, int *, double *, int *);

void dscal(int *len, double *alpha, double *x, int *incx);
double dasum(int *len, double *x, int *incx);
int idamax(int *len, double *x, int *incx);
void dsyev(char *c1,char *c2,int *n,
double *a,int *lda, double *wr,
double *wrk,int *lwrk, int *ierr);
void zdscal_(int *n,double *alpha,std::complex<double> *x,int *incx);
void dgbmv(char *trans, int *m, int *n,
int *kl, int *ku, double *alpha, double *a,
int *lda, double *x, int *incx, double *beta,
double *y, int *incy);
void dsbmv(char *uplo, int *n, int *k,
double *alpha, double *a, int *lda, double *x, int *incx,
double *beta, double *y, int *incy);
void sspev(char *vec,char *uplo,int *size,double *ap,
double *wr,double *z,int *n,double *wrk,int *ierr);
void dgesv(int *n, int *nrhs, double *a, int *lda, int *ipiv,
double *b, int *ldb, int *info);

void dvea(int*,double*,int*,double*,int*,double*,int*);
void dyax(int*,double*,double*,int*,double*,int*);
void dnaxpy(int*,int*,double*,int*,double*,int*,int*,double*,int*,int*);

#ifdef __cplusplus
}
#endif

#endif
