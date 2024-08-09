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
// MLWFTransform.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>

#include "MLWFTransform.h"
#include "D3vector.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "UnitCell.h"
#include "jade.h"
#include "blas.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MLWFTransform::MLWFTransform(const SlaterDet& sd) : sd_(sd),
 cell_(sd.basis().cell()), ctxt_(sd.context()),
 bm_(BasisMapping(sd.basis(),sd.basis().np(0),sd.basis().np(1),
 sd.basis().np(2))), maxsweep_(50), tol_(1.e-8)
{
  a_.resize(6);
  b_.resize(6);
  adiag_.resize(6);
  const int n = sd.c().n();
  const int nb = sd.c().nb();
  for ( int k = 0; k < 6; k++ )
  {
    a_[k] = new DoubleMatrix(ctxt_,n,n,nb,nb);
    b_[k] = new DoubleMatrix(ctxt_,n,n,nb,nb);
    adiag_[k].resize(n);
  }
  u_ = new DoubleMatrix(ctxt_,n,n,nb,nb);

  sdcosx_ = new SlaterDet(sd_);
  sdcosy_ = new SlaterDet(sd_);
  sdcosz_ = new SlaterDet(sd_);
  sdsinx_ = new SlaterDet(sd_);
  sdsiny_ = new SlaterDet(sd_);
  sdsinz_ = new SlaterDet(sd_);

}

////////////////////////////////////////////////////////////////////////////////
MLWFTransform::~MLWFTransform(void)
{
  for ( int k = 0; k < 6; k++ )
  {
    delete a_[k];
    delete b_[k];
  }
  delete u_;

  delete sdcosx_;
  delete sdcosy_;
  delete sdcosz_;
  delete sdsinx_;
  delete sdsiny_;
  delete sdsinz_;
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::update(void)
{
  // recompute cos and sin matrices
  const ComplexMatrix& c = sd_.c();
  ComplexMatrix& ccosx = sdcosx_->c();
  ComplexMatrix& csinx = sdsinx_->c();
  ComplexMatrix& ccosy = sdcosy_->c();
  ComplexMatrix& csiny = sdsiny_->c();
  ComplexMatrix& ccosz = sdcosz_->c();
  ComplexMatrix& csinz = sdsinz_->c();
  // proxy real matrices cr, cc, cs
  DoubleMatrix cr(c);
  DoubleMatrix ccx(ccosx);
  DoubleMatrix csx(csinx);
  DoubleMatrix ccy(ccosy);
  DoubleMatrix csy(csiny);
  DoubleMatrix ccz(ccosz);
  DoubleMatrix csz(csinz);
  vector<complex<double> > zvec(bm_.zvec_size()),
    zvec_cos(bm_.zvec_size()), zvec_sin(bm_.zvec_size()),
    ct(bm_.np012loc()), ct_cos(bm_.np012loc()), ct_sin(bm_.np012loc());

  for ( int i = 0; i < 6; i++ )
  {
    a_[i]->resize(c.n(), c.n(), c.nb(), c.nb());
    b_[i]->resize(c.n(), c.n(), c.nb(), c.nb());
    adiag_[i].resize(c.n());
  }
  u_->resize(c.n(), c.n(), c.nb(), c.nb());

  // loop over all local states
  const int np0 = bm_.np0();
  const int np1 = bm_.np1();
  const int np2 = bm_.np2();
  const int np01 = np0 * np1;
  const int np2loc = bm_.np2_loc();
  const int nvec = bm_.nvec();
  for ( int n = 0; n < c.nloc(); n++ )
  {
    const complex<double>* f = c.cvalptr(n*c.mloc());
    complex<double>* fcx = ccosx.valptr(n*c.mloc());
    complex<double>* fsx = csinx.valptr(n*c.mloc());
    complex<double>* fcy = ccosy.valptr(n*c.mloc());
    complex<double>* fsy = csiny.valptr(n*c.mloc());
    complex<double>* fcz = ccosz.valptr(n*c.mloc());
    complex<double>* fsz = csinz.valptr(n*c.mloc());

    // direction z
    // map state to array zvec_
    bm_.vector_to_zvec(&f[0],&zvec[0]);

    for ( int ivec = 0; ivec < nvec; ivec++ )
    {
      const int ibase = ivec * np2;
      compute_sincos(np2,&zvec[ibase],&zvec_cos[ibase],&zvec_sin[ibase]);
    }
    // map back zvec_cos to sdcos and zvec_sin to sdsin
    bm_.zvec_to_vector(&zvec_cos[0],&fcz[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsz[0]);

    // x direction
    // map zvec to ct
    bm_.transpose_bwd(&zvec[0],&ct[0]);

    for ( int iz = 0; iz < np2loc; iz++ )
    {
      for ( int iy = 0; iy < np1; iy++ )
      {
        const int ibase = iz * np01 + iy * np0;
        compute_sincos(np0,&ct[ibase],&ct_cos[ibase],&ct_sin[ibase]);
      }
    }
    // transpose ct_cos to zvec_cos
    bm_.transpose_fwd(&ct_cos[0],&zvec_cos[0]);
    // transpose ct_sin to zvec_sin
    bm_.transpose_fwd(&ct_sin[0],&zvec_sin[0]);

    // map back zvec_cos to sdcos and zvec_sin to sdsin
    bm_.zvec_to_vector(&zvec_cos[0],&fcx[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsx[0]);

    // y direction
    vector<complex<double> > c_tmp(np1),ccos_tmp(np1),csin_tmp(np1);
    int one = 1;
    int len = np1;
    int stride = np0;
    for ( int iz = 0; iz < np2loc; iz++ )
    {
      for ( int ix = 0; ix < np0; ix++ )
      {
        const int ibase = iz * np01 + ix;
        zcopy(&len,&ct[ibase],&stride,&c_tmp[0],&one);
        compute_sincos(np1,&c_tmp[0],&ccos_tmp[0],&csin_tmp[0]);
        zcopy(&len,&ccos_tmp[0],&one,&ct_cos[ibase],&stride);
        zcopy(&len,&csin_tmp[0],&one,&ct_sin[ibase],&stride);
      }
    }
    // transpose ct_cos to zvec_cos
    bm_.transpose_fwd(&ct_cos[0],&zvec_cos[0]);
    // transpose ct_sin to zvec_sin
    bm_.transpose_fwd(&ct_sin[0],&zvec_sin[0]);

    // map back zvec_cos and zvec_sin
    bm_.zvec_to_vector(&zvec_cos[0],&fcy[0]);
    bm_.zvec_to_vector(&zvec_sin[0],&fsy[0]);
  }

  // dot products b_[0] = <cos x>, b_[1] = <sin x>
  b_[0]->gemm('t','n',2.0,cr,ccx,0.0);
  b_[0]->ger(-1.0,cr,0,ccx,0);
  b_[1]->gemm('t','n',2.0,cr,csx,0.0);
  b_[1]->ger(-1.0,cr,0,csx,0);

  // dot products b_[2] = <cos y>, b_[3] = <sin y>
  b_[2]->gemm('t','n',2.0,cr,ccy,0.0);
  b_[2]->ger(-1.0,cr,0,ccy,0);
  b_[3]->gemm('t','n',2.0,cr,csy,0.0);
  b_[3]->ger(-1.0,cr,0,csy,0);

  // dot products b_[4] = <cos z>, b_[5] = <sin z>
  b_[4]->gemm('t','n',2.0,cr,ccz,0.0);
  b_[4]->ger(-1.0,cr,0,ccz,0);
  b_[5]->gemm('t','n',2.0,cr,csz,0.0);
  b_[5]->ger(-1.0,cr,0,csz,0);

  // a_ = A * b_
  const double *amat = cell_.amat();
  // a_[0]
  a_[0]->clear();
  a_[0]->axpy(amat[0], *b_[0]);
  a_[0]->axpy(amat[3], *b_[2]);
  a_[0]->axpy(amat[6], *b_[4]);

  // a_[1]
  a_[1]->clear();
  a_[1]->axpy(amat[0], *b_[1]);
  a_[1]->axpy(amat[3], *b_[3]);
  a_[1]->axpy(amat[6], *b_[5]);

  // a_[2]
  a_[2]->clear();
  a_[2]->axpy(amat[1], *b_[0]);
  a_[2]->axpy(amat[4], *b_[2]);
  a_[2]->axpy(amat[7], *b_[4]);

  // a_[3]
  a_[3]->clear();
  a_[3]->axpy(amat[1], *b_[1]);
  a_[3]->axpy(amat[4], *b_[3]);
  a_[3]->axpy(amat[7], *b_[5]);

  // a_[4]
  a_[4]->clear();
  a_[4]->axpy(amat[2], *b_[0]);
  a_[4]->axpy(amat[5], *b_[2]);
  a_[4]->axpy(amat[8], *b_[4]);

  // a_[5]
  a_[5]->clear();
  a_[5]->axpy(amat[2], *b_[1]);
  a_[5]->axpy(amat[5], *b_[3]);
  a_[5]->axpy(amat[8], *b_[5]);
}
////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::compute_transform(void)
{
  jade(maxsweep_,tol_,a_,*u_,adiag_);
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::compute_sincos(const int n, const complex<double>* f,
  complex<double>* fc, complex<double>* fs)
{
  // fc[i] =     0.5 * ( f[i-1] + f[i+1] )
  // fs[i] = (0.5/i) * ( f[i-1] - f[i+1] )

  // i = 0
  complex<double> zp = f[n-1];
  complex<double> zm = f[1];
  fc[0] = 0.5 * ( zp + zm );
  complex<double> zdiff = zp - zm;
  fs[0] = 0.5 * complex<double>(imag(zdiff),-real(zdiff));
  for ( int i = 1; i < n-1; i++ )
  {
    const complex<double> zzp = f[i-1];
    const complex<double> zzm = f[i+1];
    fc[i] = 0.5 * ( zzp + zzm );
    const complex<double> zzdiff = zzp - zzm;
    fs[i] = 0.5 * complex<double>(imag(zzdiff),-real(zzdiff));
  }
  // i = n-1
  zp = f[n-2];
  zm = f[0];
  fc[n-1] = 0.5 * ( zp + zm );
  zdiff = zp - zm;
  fs[n-1] = 0.5 * complex<double>(imag(zdiff),-real(zdiff));
}

////////////////////////////////////////////////////////////////////////////////
D3vector MLWFTransform::center(int i) const
{
  assert(i>=0 && i<sd_.nst());
  // c,s = B^T * adiag
  const double *bmat = cell_.bmat();
  const double c0 = bmat[0] * adiag_[0][i] +
                    bmat[1] * adiag_[2][i] +
                    bmat[2] * adiag_[4][i];
  const double s0 = bmat[0] * adiag_[1][i] +
                    bmat[1] * adiag_[3][i] +
                    bmat[2] * adiag_[5][i];

  const double c1 = bmat[3] * adiag_[0][i] +
                    bmat[4] * adiag_[2][i] +
                    bmat[5] * adiag_[4][i];
  const double s1 = bmat[3] * adiag_[1][i] +
                    bmat[4] * adiag_[3][i] +
                    bmat[5] * adiag_[5][i];

  const double c2 = bmat[6] * adiag_[0][i] +
                    bmat[7] * adiag_[2][i] +
                    bmat[8] * adiag_[4][i];
  const double s2 = bmat[6] * adiag_[1][i] +
                    bmat[7] * adiag_[3][i] +
                    bmat[8] * adiag_[5][i];

  const double itwopi = 1.0 / ( 2.0 * M_PI );
  const double t0 = itwopi * atan2(s0,c0);
  const double t1 = itwopi * atan2(s1,c1);
  const double t2 = itwopi * atan2(s2,c2);
  const double *amat = cell_.amat();
  const double x = t0 * amat[0] + t1 * amat[3] + t2 * amat[6];
  const double y = t0 * amat[1] + t1 * amat[4] + t2 * amat[7];
  const double z = t0 * amat[2] + t1 * amat[5] + t2 * amat[8];
  return D3vector(x,y,z);
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(int i, int j) const
{
  // squared spread of state i in the direction of
  // the reciprocal lattice vector b_j
  assert(i>=0 && i<sd_.nst());
  assert(j>=0 && j<3);

  const double itwopi = 1.0 / ( 2.0 * M_PI );
  const double *bmat = cell_.bmat();
  // c,s = B^T * adiag / ( 2 * pi )
  const double c = itwopi * ( bmat[3*j+0] * adiag_[0][i] +
                              bmat[3*j+1] * adiag_[2][i] +
                              bmat[3*j+2] * adiag_[4][i] );
  const double s = itwopi * ( bmat[3*j+0] * adiag_[1][i] +
                              bmat[3*j+1] * adiag_[3][i] +
                              bmat[3*j+2] * adiag_[5][i] );

  const double fac = 1.0 / length(cell_.b(j));
  return fac*fac * ( 1.0 - c*c - s*s );
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(int i) const
{
  assert(i>=0 & i<sd_.nst());
  return spread2(i,0) + spread2(i,1) + spread2(i,2);
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread(int i) const
{
  return sqrt(spread2(i));
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread2(void) const
{
  double sum = 0.0;
  for ( int i = 0; i < sd_.nst(); i++ )
    sum += spread2(i);
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double MLWFTransform::spread(void) const
{
  return sqrt(spread2());
}

////////////////////////////////////////////////////////////////////////////////
D3vector MLWFTransform::dipole(void) const
{
  // total electronic dipole
  D3vector sum(0.0,0.0,0.0);
  for ( int i = 0; i < sd_.nst(); i++ )
    sum -= sd_.occ(i) * center(i);
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::apply_transform(SlaterDet& sd)
{
  // proxy double matrix c
  DoubleMatrix c(sd.c());
  DoubleMatrix cp(c);
  c.gemm('n','n',1.0,cp,*u_,0.0);
}

////////////////////////////////////////////////////////////////////////////////
void MLWFTransform::print(ostream& os) const
{
  for ( int i = 0; i < sd_.nst(); i++ )
  {
    D3vector ctr = center(i);

    os.setf(ios::fixed, ios::floatfield);
    os.setf(ios::right, ios::adjustfield);
    os << "   <mlwf>\n"
       << "     <center>  " << setprecision(8)
       << setw(14) << ctr.x
       << setw(14) << ctr.y
       << setw(14) << ctr.z
       << " </center>\n"
       << "     <spread2> "
       << setw(14) << spread2(i,0)
       << setw(14) << spread2(i,1)
       << setw(14) << spread2(i,2)
       << " </spread2>\n"
       << "   </mlwf>"
       << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, MLWFTransform& mlwft)
{
  mlwft.print(os);
  return os;
}
