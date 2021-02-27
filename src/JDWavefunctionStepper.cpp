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
// JDWavefunctionStepper.cpp
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
  Preconditioner& prec, EnergyFunctional& ef, TimerMap& tmap) :
  WavefunctionStepper(wf,tmap), prec_(prec), wft_(wf), dwft_(wf), ef_(ef)
{
  tmap_["jd_residual"].reset();
  tmap_["jd_compute_z"].reset();
  tmap_["jd_hz"].reset();
  tmap_["jd_blocks"].reset();
  tmap_["jd_gemm"].reset();
  tmap_["jd_getsub"].reset();
  tmap_["jd_syev"].reset();
  tmap_["jd_heev"].reset();
}

////////////////////////////////////////////////////////////////////////////////
JDWavefunctionStepper::~JDWavefunctionStepper(void)
{}

////////////////////////////////////////////////////////////////////////////////
void JDWavefunctionStepper::update(Wavefunction& dwf)
{
  // save copy of wf in wft_ and dwf in dwft_
  wft_ = wf_;
  dwft_ = dwf;
  tmap_["jd_residual"].start();
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      SlaterDet* sd = wf_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd = dwf.sd(isp_loc,ikp_loc);
      if ( sd->basis().real() )
      {
        // compute A = Y^T H Y  and descent direction HY - YA
        // proxy real matrices c, cp
        DoubleMatrix c_proxy(sd->c());
        DoubleMatrix cp_proxy(dsd->c());
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
        ComplexMatrix& c = sd->c();
        ComplexMatrix& cp = dsd->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());

        // (Y,HY)
        a.gemm('c','n',1.0,c,cp,0.0);

        // cp = cp - c * a
        cp.gemm('n','n',-1.0,c,a,1.0);
        // dwf.sd->c() now contains the descent direction (HV-VA)
      }
    } // ikp
  } // isp_loc
  tmap_["jd_residual"].stop();

  // dwf.sd->c() now contains the descent direction (HV-VA)

  // update preconditioner
  prec_.update(wf_);

  tmap_["jd_compute_z"].start();
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      if ( wf_.sd(isp_loc,ikp_loc)->basis().real() )
      {
        // Apply preconditioner K and store -K(HV-VA) in wf

        double* c = (double*) wf_.sd(isp_loc,ikp_loc)->c().valptr();
        double* dc = (double*) dwf.sd(isp_loc,ikp_loc)->c().valptr();
        DoubleMatrix c_proxy(wf_.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());
        DoubleMatrix c_proxy_t(wft_.sd(isp_loc,ikp_loc)->c());
        const int mloc = wf_.sd(isp_loc,ikp_loc)->c().mloc();
        const int ngwl = wf_.sd(isp_loc,ikp_loc)->basis().localsize();
        const int nloc = wf_.sd(isp_loc,ikp_loc)->c().nloc();

        for ( int n = 0; n < nloc; n++ )
        {
          // note: double mloc length for complex<double> indices
          double* dcn = &dc[2*mloc*n];
          double* cn  =  &c[2*mloc*n];

          for ( int i = 0; i < ngwl; i++ )
          {
            const double fac = prec_.diag(isp_loc,ikp_loc,n,i);
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
        wf_.sd(isp_loc,ikp_loc)->gram();

        // wf now contains Z, orthonormal
      } // if real
      else
      {
        // Apply preconditioner K and store -K(HV-VA) in wf
        ComplexMatrix& c = wf_.sd(isp_loc,ikp_loc)->c();
        complex<double>* cv = c.valptr();
        ComplexMatrix& cp = dwf.sd(isp_loc,ikp_loc)->c();
        complex<double>* cpv = cp.valptr();
        ComplexMatrix& ct = wft_.sd(isp_loc,ikp_loc)->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        const int ngwl = wf_.sd(isp_loc,ikp_loc)->basis().localsize();
        const int mloc = c.mloc();
        const int nloc = c.nloc();

        for ( int n = 0; n < nloc; n++ )
        {
          complex<double>* cpn = &cpv[mloc*n];
          complex<double>* cn  =  &cv[mloc*n];

          for ( int i = 0; i < ngwl; i++ )
          {
            const double fac = prec_.diag(isp_loc,ikp_loc,n,i);
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
        wf_.sd(isp_loc,ikp_loc)->gram();

        // wf now contains Z, orthonormal
      }
    } // ikp_loc
  } // isp_loc
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
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      if ( wf_.sd(isp_loc,ikp_loc)->basis().real() )
      {
        // wf now contains Z
        // dwf now contains HZ
        // compute blocks (Y,HY), (Y,HZ), (Z,HZ)

        // proxy real matrices c, cp
        DoubleMatrix c_proxy(wf_.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix c_proxy_t(wft_.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix cp_proxy(dwf.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix cp_proxy_t(dwft_.sd(isp_loc,ikp_loc)->c());

        DoubleMatrix a(c_proxy.context(),c_proxy.n(),c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());
        DoubleMatrix h(c_proxy.context(),2*c_proxy.n(),2*c_proxy.n(),
                       c_proxy.nb(),c_proxy.nb());

        // (Y,HY)
        // factor 2.0 in next line: G and -G
        tmap_["jd_gemm"].start();
        a.gemm('t','n',2.0,c_proxy_t,cp_proxy_t,0.0);
        tmap_["jd_gemm"].stop();
        // rank-1 update correction
        a.ger(-1.0,c_proxy_t,0,cp_proxy_t,0);
        // a contains (Y,HY), copy to h11 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,0,0);
        tmap_["jd_getsub"].stop();

        // (Z,HY)
        // factor 2.0 in next line: G and -G
        tmap_["jd_gemm"].start();
        a.gemm('t','n',2.0,c_proxy,cp_proxy_t,0.0);
        tmap_["jd_gemm"].stop();
        // rank-1 update correction
        a.ger(-1.0,c_proxy,0,cp_proxy_t,0);
        // a contains (Z,HY), copy to h21 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,a.m(),0);
        tmap_["jd_getsub"].stop();

        // (Z,HZ)
        // factor 2.0 in next line: G and -G
        tmap_["jd_gemm"].start();
        a.gemm('t','n',2.0,c_proxy,cp_proxy,0.0);
        tmap_["jd_gemm"].stop();
        // rank-1 update correction
        a.ger(-1.0,c_proxy,0,cp_proxy,0);
        // a contains (Z,HZ), copy to h22 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,a.m(),a.n());
        tmap_["jd_getsub"].stop();

        // diagonalize h
        // Note: we only need the first n eigenvectors of the (2n x 2n) matrix
        valarray<double> w(h.m());
        // q is (2n,2n)
        DoubleMatrix q(h.context(),h.n(),h.n(),h.nb(),h.nb());
        tmap_["jd_syev"].start();
        h.syevd('l',w,q);
        tmap_["jd_syev"].stop();

        // compute the first n eigenvectors and store in wf
        // Y = Z Q21 (store result in dwf)
        // get Q21 in a
        tmap_["jd_getsub"].start();
        a.getsub(q,a.n(),a.n(),a.n(),0);
        tmap_["jd_getsub"].stop();
        tmap_["jd_gemm"].start();
        cp_proxy.gemm('n','n',1.0,c_proxy,a,0.0);
        tmap_["jd_gemm"].stop();

        // Y = Y Q11 (store result in wf)
        // get Q11 in a
        tmap_["jd_getsub"].start();
        a.getsub(q,a.n(),a.n(),0,0);
        tmap_["jd_getsub"].stop();
        tmap_["jd_gemm"].start();
        c_proxy.gemm('n','n',1.0,c_proxy_t,a,0.0);
        tmap_["jd_gemm"].stop();

        // add two contributions
        c_proxy += cp_proxy;

        // wf now contains the corrected eigenvectors
      } // if real
      else
      {
        // wf now contains Z
        // dwf now contains HZ
        // compute blocks (Y,HY), (Y,HZ), (Z,HZ)

        ComplexMatrix& c = wf_.sd(isp_loc,ikp_loc)->c();
        ComplexMatrix& ct = wft_.sd(isp_loc,ikp_loc)->c();
        ComplexMatrix& cp = dwf.sd(isp_loc,ikp_loc)->c();
        ComplexMatrix& cpt = dwft_.sd(isp_loc,ikp_loc)->c();

        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        ComplexMatrix h(c.context(),2*c.n(),2*c.n(),c.nb(),c.nb());

        // (Y,HY)
        // factor 2.0 in next line: G and -G
        tmap_["jd_gemm"].start();
        a.gemm('c','n',1.0,ct,cpt,0.0);
        tmap_["jd_gemm"].stop();
        // a contains (Y,HY), copy to h11 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,0,0);
        tmap_["jd_getsub"].stop();

        // (Z,HY)
        tmap_["jd_gemm"].start();
        a.gemm('c','n',1.0,c,cpt,0.0);
        tmap_["jd_gemm"].stop();
        // a contains (Z,HY), copy to h21 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,a.m(),0);
        tmap_["jd_getsub"].stop();

        // (Z,HZ)
        tmap_["jd_gemm"].start();
        a.gemm('c','n',1.0,c,cp,0.0);
        tmap_["jd_gemm"].stop();
        // a contains (Z,HZ), copy to h22 block
        tmap_["jd_getsub"].start();
        h.getsub(a,a.m(),a.n(),0,0,a.m(),a.n());
        tmap_["jd_getsub"].stop();

        // diagonalize h
        // Note: we only need the first n eigenvectors of the (2n x 2n) matrix
        valarray<double> w(h.m());
        // q is (2n,2n)
        ComplexMatrix q(h.context(),h.n(),h.n(),h.nb(),h.nb());
        tmap_["jd_heev"].start();
        h.heevd('l',w,q);
        tmap_["jd_heev"].stop();

        // compute the first n eigenvectors and store in wf
        // Y = Z Q21 (store result in dwf)
        // get Q21 in a
        tmap_["jd_getsub"].start();
        a.getsub(q,a.n(),a.n(),a.n(),0);
        tmap_["jd_getsub"].stop();
        tmap_["jd_gemm"].start();
        cp.gemm('n','n',1.0,c,a,0.0);
        tmap_["jd_gemm"].stop();

        // Y = Y Q11 (store result in wf)
        // get Q11 in a
        tmap_["jd_getsub"].start();
        a.getsub(q,a.n(),a.n(),0,0);
        tmap_["jd_getsub"].stop();
        tmap_["jd_gemm"].start();
        c.gemm('n','n',1.0,ct,a,0.0);
        tmap_["jd_gemm"].stop();

        // add two contributions
        c += cp;

        // wf now contains the corrected eigenvectors
      } // if real
    } // ikp_loc
  } // isp_loc
  tmap_["jd_blocks"].stop();
}
