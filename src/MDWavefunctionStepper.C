////////////////////////////////////////////////////////////////////////////////
//
// MDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDWavefunctionStepper.C,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#include "MDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
MDWavefunctionStepper::MDWavefunctionStepper(Sample& s, TimerMap& tmap) : 
  s_(s), wf_(s.wf), tmap_(tmap)
{
  dt_ = s_.ctrl.dt;
  const double emass = s_.ctrl.emass;
  dt2bye_ = (emass == 0.0) ? 0.5 / wf_.ecut() : dt_*dt_/emass;           
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::update(Wavefunction& dwf)
{
  // Verlet update of wf using force dwf and wfm stored in *wfv
  Wavefunction* const wfv = s_.wfv;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          tmap_["update_psi"].start();
          // Verlet update of wf
          // cp = c + (c - cm) + dt2/m * hpsi
          // c += c - cm + dt2bye * hpsi
          // cm = c
          SlaterDet* sd = wf_.sd(ispin,ikp);
          double* cptr = (double*) sd->c().valptr();
          double* cptrm = (double*) wfv->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
            (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cm = &cptrm[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cm[2*i];
              const double cmtmp1 = cm[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];
              const double cptmp = 2.0*ctmp -  cmtmp -  dt2bye_ * dctmp;
              const double cptmp1= 2.0*ctmp1 - cmtmp1 - dt2bye_ * dctmp1;
 
              c[2*i]    = cptmp;
              c[2*i+1]  = cptmp1;
              cm[2*i]   = ctmp;
              cm[2*i+1] = ctmp1;
            }
          }
          tmap_["update_psi"].stop();
          
          tmap_["riccati"].start();
          assert(wfv!=0);
          wf_.sd(ispin,ikp)->riccati(*wfv->sd(ispin,ikp));
          tmap_["riccati"].stop();
        }
      }
    }
  }
  ekin_em_ = ekin_ep_;
  ekin_ep_ = ekin_eh();
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::stoermer_start(Wavefunction& dwf)
{
  // First step MD step using Stoermer rule
  
  Wavefunction* const wfv = s_.wfv;
  assert(wfv!=0);
  // First iteration of Stoermer's rule
  // cp = c + dt * v - 0.5 * dt2/m * hpsi
  // replace wfv by wfm
  const double m_e = dt_ * dt_ / dt2bye_;
  const double half_dt2bye = 0.5 * dt2bye_;
  double ekin_e = 0.0;
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          tmap_["update_psi"].start();
          SlaterDet* sd = wf_.sd(ispin,ikp);
          double* cptr = (double*) sd->c().valptr();
          double* cptrv = (double*) wfv->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
            (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const vector<double>& occ = sd->occ();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          const bool onrow0 = ( wf_.context().myrow() == 0 );
          for ( int n = 0; n < nloc; n++ )
          {
            const int nglobal = sd->c().j(0,n);
            const double occn = occ[nglobal];
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cv = &cptrv[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            double tmpsum = 0.0;
            if ( onrow0 )
            {
              const double cvtmp = cv[0];
              tmpsum -= 0.5 * (cvtmp*cvtmp);
              cv[1] = 0.0;
            }
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cvtmp = cv[2*i];
              const double cvtmp1 = cv[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];

              tmpsum += (cvtmp*cvtmp + cvtmp1*cvtmp1);
              c[2*i]    = ctmp  + dt_ * cvtmp  - half_dt2bye * dctmp;
              c[2*i+1]  = ctmp1 + dt_ * cvtmp1 - half_dt2bye * dctmp1;
              cv[2*i]   = ctmp;
              cv[2*i+1] = ctmp1;
              // Note: 2 in next line from G,-G
            }
            ekin_e += 2.0 * m_e * occn * tmpsum;
          }
          tmap_["update_psi"].stop();
          
          tmap_["riccati"].start();
          assert(wfv!=0);
          wf_.sd(ispin,ikp)->riccati(*wfv->sd(ispin,ikp));
          tmap_["riccati"].stop();
        }
      }
    }
  }
  wf_.context().dsum(1,1,&ekin_e,1);
  ekin_em_ = ekin_ep_ = ekin_e;
  // Note: *wfv now contains wf(t-dt)
}

////////////////////////////////////////////////////////////////////////////////
void MDWavefunctionStepper::stoermer_end(Wavefunction& dwf)
{
  // Last step of Stoermer rule
  // Compute wfv = (wf - wfm)/dt + 0.5*dtbye*dwf
  
  Wavefunction * const wfv = s_.wfv;
  assert(wfv!=0);
  assert(dt_!=0.0);
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          // Last iteration of Stoermer's rule
          // compute final velocity wfv
          // v = ( c - cm ) / dt - 0.5 * dt/m * hpsi
 
          // Note: At this point, *wfv contains wf(t-dt)
 
          // hpsi must be orthogonal to subspace spanned by c
          // compute descent direction H psi - psi (psi^T H psi)
 
          SlaterDet* sd = wf_.sd(ispin,ikp);
          if ( sd->basis().real() )
          {
            // proxy real matrices c, cp
            DoubleMatrix c(sd->c());
            DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
            DoubleMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
            // factor 2.0 in next line: G and -G
            a.gemm('t','n',2.0,c,cp,0.0);
            // rank-1 update correction
            a.ger(-1.0,c,0,cp,0);
 
            // cp = cp - c * a
            cp.gemm('n','n',-1.0,c,a,1.0);
          }
          else
          {
            // not implemented in the complex case
            assert(false);
          }
 
          const double dt_inv = 1.0/dt_;
          const double half_dtbye = 0.5 * dt2bye_ / dt_;
          double* cptr = (double*) sd->c().valptr();
          double* cptrv = (double*) wfv->sd(ispin,ikp)->c().valptr();
          const double* dcptr =
            (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cv = &cptrv[2*mloc*n];
            const double* dc = &dcptr[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cv[2*i];
              const double cmtmp1 = cv[2*i+1];
              const double dctmp = dc[2*i];
              const double dctmp1 = dc[2*i+1];
 
              cv[2*i]   = ( ctmp  - cmtmp  ) * dt_inv
                          - half_dtbye * dctmp;
              cv[2*i+1] = ( ctmp1 - cmtmp1 ) * dt_inv
                          - half_dtbye * dctmp1;
            }
          }
          // Note: *wfv now contains the wavefunction velocity
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double MDWavefunctionStepper::ekin_eh(void)
{
  // compute ekin at time t - 0.5*dt using wf and wfm
  tmap_["ekin_e"].start();
  double ekin_e = 0.0;
  Wavefunction* const wfv = s_.wfv; // assume that wfv contains wfm
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          SlaterDet* sd = wf_.sd(ispin,ikp);
          double* cptr = (double*) sd->c().valptr();
          double* cptrm = (double*) wfv->sd(ispin,ikp)->c().valptr();
          const int mloc = sd->c().mloc();
          const int nloc = sd->c().nloc();
          // compute electronic kinetic energy at time t-1/2
          const bool onrow0 = ( wf_.context().myrow() == 0 );
          const vector<double>& occ = sd->occ();
          for ( int n = 0; n < nloc; n++ )
          {
            const int nglobal = sd->c().j(0,n);
            const double occn = occ[nglobal];
            // note: double mloc length for complex<double> indices
            double* c = &cptr[2*mloc*n];
            double* cm = &cptrm[2*mloc*n];
            double tmpsum = 0.0;
            if ( onrow0 )
            {
              // correct for double counting of G=0 element
              // i=0 coefficient is real, count only real part
              const double ctmp = c[0];
              const double cmtmp = cm[0];
              tmpsum -= 0.5 * (ctmp - cmtmp)*(ctmp - cmtmp);
            }
            for ( int i = 0; i < mloc; i++ )
            {
              const double ctmp = c[2*i];
              const double ctmp1 = c[2*i+1];
              const double cmtmp = cm[2*i];
              const double cmtmp1 = cm[2*i+1];
 
              tmpsum += (ctmp -cmtmp )*(ctmp -cmtmp ) +
                        (ctmp1-cmtmp1)*(ctmp1-cmtmp1);
            }
            // Note: 2 in next line: from (G,-G)
            ekin_e += ( 2.0 * occn / dt2bye_ ) * tmpsum;
          }
        }
      }
    }
  }
  wf_.context().dsum(1,1,&ekin_e,1);
  tmap_["ekin_e"].stop();
  return ekin_e;
}
