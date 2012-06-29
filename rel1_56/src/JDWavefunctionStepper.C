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
// JDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "JDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "EnergyFunctional.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
JDWavefunctionStepper::JDWavefunctionStepper(Wavefunction& wf,
  double ecutprec, EnergyFunctional& ef, TimerMap& tmap) :
  WavefunctionStepper(wf,tmap), prec_(0), wft_(wf), dwft_(wf), ef_(ef)
{
  prec_ = new Preconditioner(wf,ecutprec);
}

////////////////////////////////////////////////////////////////////////////////
JDWavefunctionStepper::~JDWavefunctionStepper(void)
{
  delete prec_;
}

////////////////////////////////////////////////////////////////////////////////
void JDWavefunctionStepper::update(Wavefunction& dwf)
{
  // save copy of wf in wft_ and dwf in dwft_
  wft_ = wf_;
  dwft_ = dwf;
  tmap_["jd_residual"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp)->basis().real() )
      {
        // compute A = Y^T H Y  and descent direction HY - YA
        // proxy real matrices c, cp
        DoubleMatrix c_proxy(wf_.sd(ispin,ikp)->c());
        DoubleMatrix cp_proxy(dwf.sd(ispin,ikp)->c());
        DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());

        // factor 2.0 in next line: G and -G
        a.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
        // rank-1 update correction
        a.ger(-1.0,c_proxy,0,cp_proxy,0);

        // cp = cp - c * a
        cp_proxy.gemm('n','n',-1.0,c_proxy,a,1.0);
      }
      else // real
      {
        // compute A = Y^T H Y  and descent direction HY - YA
        ComplexMatrix& c = wf_.sd(ispin,ikp)->c();
        ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());

        // (Y,HY)
        a.gemm('c','n',1.0,c,cp,0.0);

        // cp = cp - c * a
        cp.gemm('n','n',-1.0,c,a,1.0);
        // dwf.sd->c() now contains the descent direction (HV-VA)
      }
    } // ikp
  } // ispin
  tmap_["jd_residual"].stop();

  // dwf.sd->c() now contains the descent direction (HV-VA)
  prec_->update(wf_);

  tmap_["jd_compute_z"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp)->basis().real() )
      {
        // Apply preconditioner K and store -K(HV-VA) in wf

        double* c = (double*) wf_.sd(ispin,ikp)->c().valptr();
        double* dc = (double*) dwf.sd(ispin,ikp)->c().valptr();
        DoubleMatrix c_proxy(wf_.sd(ispin,ikp)->c());
        DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());
        DoubleMatrix c_proxy_t(wft_.sd(ispin,ikp)->c());
        const int mloc = wf_.sd(ispin,ikp)->c().mloc();
        const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
        const int nloc = wf_.sd(ispin,ikp)->c().nloc();

        for ( int n = 0; n < nloc; n++ )
        {
          // note: double mloc length for complex<double> indices
          double* dcn = &dc[2*mloc*n];
          double* cn  =  &c[2*mloc*n];

          for ( int i = 0; i < ngwl; i++ )
          {
            const double fac = prec_->diag(ispin,ikp,n,i);
            const double f0 = -fac * dcn[2*i];
            const double f1 = -fac * dcn[2*i+1];
            cn[2*i] = f0;
            cn[2*i+1] = f1;
          }
        }
        // wf_ now contains the preconditioned descent
        // direction Z = -K(HY-YA)

        // orthogonalize Z to Y
        // Z = Z - YY^T Z
        // A = Y^T * Z
        // factor 2.0 in next line: G and -G
        a.gemm('t','n',2.0,c_proxy_t,c_proxy,0.0);
        // rank-1 update correction
        a.ger(-1.0,c_proxy_t,0,c_proxy,0);

        // Z = Z - Y * A
        c_proxy.gemm('n','n',-1.0,c_proxy_t,a,1.0);

        // orthogonalize Z: gram(Z)
        wf_.sd(ispin,ikp)->gram();

        // wf now contains Z, orthonormal
      } // if real
      else
      {
        // Apply preconditioner K and store -K(HV-VA) in wf
        ComplexMatrix& c = wf_.sd(ispin,ikp)->c();
        complex<double>* cv = c.valptr();
        ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
        complex<double>* cpv = cp.valptr();
        ComplexMatrix& ct = wft_.sd(ispin,ikp)->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
        const int mloc = c.mloc();
        const int nloc = c.nloc();

        for ( int n = 0; n < nloc; n++ )
        {
          complex<double>* cpn = &cpv[mloc*n];
          complex<double>* cn  =  &cv[mloc*n];

          for ( int i = 0; i < ngwl; i++ )
          {
            const double fac = prec_->diag(ispin,ikp,n,i);
            cn[i] = -fac * cpn[i];
          }
        }
        // wf_ now contains the preconditioned descent
        // direction Z = -K(HY-YA)

        // orthogonalize Z to Y
        // Z = Z - YY^T Z
        // A = Y^T * Z
        a.gemm('c','n',1.0,ct,c,0.0);

        // Z = Z - Y * A
        c.gemm('n','n',-1.0,ct,a,1.0);

        // orthogonalize Z: gram(Z)
        wf_.sd(ispin,ikp)->gram();

        // wf now contains Z, orthonormal
      }
    } // ikp
  } // ispin
  tmap_["jd_compute_z"].stop();

  tmap_["jd_hz"].start();
  // compute HZ
  const bool compute_hpsi = true;
  const bool compute_forces = false;
  const bool compute_stress = false;
  vector<vector<double> > fion;
  valarray<double> sigma;
  ef_.energy(compute_hpsi,dwf,compute_forces,fion,
             compute_stress,sigma);
  tmap_["jd_hz"].stop();

  tmap_["jd_blocks"].start();
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp)->basis().real() )
      {
        // wf now contains Z
        // dwf now contains HZ
        // compute blocks (Y,HY), (Y,HZ), (Z,HZ)

        // proxy real matrices c, cp
        DoubleMatrix c_proxy(wf_.sd(ispin,ikp)->c());
        DoubleMatrix c_proxy_t(wft_.sd(ispin,ikp)->c());
        DoubleMatrix cp_proxy(dwf.sd(ispin,ikp)->c());
        DoubleMatrix cp_proxy_t(dwft_.sd(ispin,ikp)->c());

        DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());
        DoubleMatrix h(c_proxy.context(),2*c_proxy.n(),2*c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());

        // (Y,HY)
        // factor 2.0 in next line: G and -G
        a.gemm('t','n',2.0,c_proxy_t,cp_proxy_t,0.0);
        // rank-1 update correction
        a.ger(-1.0,c_proxy_t,0,cp_proxy_t,0);
        // a contains (Y,HY), copy to h11 block
        h.getsub(a,a.m(),a.n(),0,0,0,0);

        // (Z,HY)
        // factor 2.0 in next line: G and -G
        a.gemm('t','n',2.0,c_proxy,cp_proxy_t,0.0);
        // rank-1 update correction
        a.ger(-1.0,c_proxy,0,cp_proxy_t,0);
        // a contains (Z,HY), copy to h21 block
        h.getsub(a,a.m(),a.n(),0,0,a.m(),0);

        // (Z,HZ)
        // factor 2.0 in next line: G and -G
        a.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
        // rank-1 update correction
        a.ger(-1.0,c_proxy,0,cp_proxy,0);
        // a contains (Z,HZ), copy to h22 block
        h.getsub(a,a.m(),a.n(),0,0,a.m(),a.n());

        // diagonalize h
        // Note: we only need the first n eigenvectors of the (2n x 2n) matrix
        valarray<double> w(h.m());
        // q is (2n,2n)
        DoubleMatrix q(h.context(),h.n(),h.n(),h.nb(),h.nb());
        h.syev('l',w,q);

        // compute the first n eigenvectors and store in wf
        // Y = Z Q21 (store result in dwf)
        // get Q21 in a
        a.getsub(q,a.n(),a.n(),a.n(),0);
        cp_proxy.gemm('n','n',1.0,c_proxy,a,0.0);

        // Y = Y Q11 (store result in wf)
        // get Q11 in a
        a.getsub(q,a.n(),a.n(),0,0);
        c_proxy.gemm('n','n',1.0,c_proxy_t,a,0.0);

        // add two contributions
        c_proxy += cp_proxy;

        // wf now contains the corrected eigenvectors
      } // if real
      else
      {
        // wf now contains Z
        // dwf now contains HZ
        // compute blocks (Y,HY), (Y,HZ), (Z,HZ)

        ComplexMatrix& c = wf_.sd(ispin,ikp)->c();
        ComplexMatrix& ct = wft_.sd(ispin,ikp)->c();
        ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
        ComplexMatrix& cpt = dwft_.sd(ispin,ikp)->c();

        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        ComplexMatrix h(c.context(),2*c.n(),2*c.n(),c.nb(),c.nb());

        // (Y,HY)
        // factor 2.0 in next line: G and -G
        a.gemm('c','n',1.0,ct,cpt,0.0);
        // a contains (Y,HY), copy to h11 block
        h.getsub(a,a.m(),a.n(),0,0,0,0);

        // (Z,HY)
        a.gemm('c','n',1.0,c,cpt,0.0);
        // a contains (Z,HY), copy to h21 block
        h.getsub(a,a.m(),a.n(),0,0,a.m(),0);

        // (Z,HZ)
        a.gemm('c','n',1.0,c,cp,0.0);
        // a contains (Z,HZ), copy to h22 block
        h.getsub(a,a.m(),a.n(),0,0,a.m(),a.n());

        // diagonalize h
        // Note: we only need the first n eigenvectors of the (2n x 2n) matrix
        valarray<double> w(h.m());
        // q is (2n,2n)
        ComplexMatrix q(h.context(),h.n(),h.n(),h.nb(),h.nb());
        h.heev('l',w,q);

        // compute the first n eigenvectors and store in wf
        // Y = Z Q21 (store result in dwf)
        // get Q21 in a
        a.getsub(q,a.n(),a.n(),a.n(),0);
        cp.gemm('n','n',1.0,c,a,0.0);

        // Y = Y Q11 (store result in wf)
        // get Q11 in a
        a.getsub(q,a.n(),a.n(),0,0);
        c.gemm('n','n',1.0,ct,a,0.0);

        // add two contributions
        c += cp;

        // wf now contains the corrected eigenvectors
      } // if real
    } // ikp
  } // ispin
  tmap_["jd_blocks"].stop();
}
