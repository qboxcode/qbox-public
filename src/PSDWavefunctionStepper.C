////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.C,v 1.2 2003-12-19 00:33:32 fgygi Exp $

#include "PSDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDWavefunctionStepper::PSDWavefunctionStepper(Sample& s, TimerMap& tmap) : 
  s_(s), wf_(s.wf), tmap_(tmap)
{
  dt_ = s_.ctrl.dt;
  const double emass = s_.ctrl.emass;
  dt2bye_ = (emass == 0.0) ? 0.5 / wf_.ecut() : dt_*dt_/emass;           
}

////////////////////////////////////////////////////////////////////////////////
void PSDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          // compute A = V^T H V  and descent direction HV - VA
 
          if ( wf_.sd(ispin,ikp)->basis().real() )
          {
            // proxy real matrices c, cp
            DoubleMatrix c(wf_.sd(ispin,ikp)->c());
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
 
          // dwf.sd->c() now contains the descent direction (HV-VA)
 
          const double g2i_prec = s_.ctrl.ecutprec > 0.0 ? 
                                  0.5 / s_.ctrl.ecutprec :
                                  0.5 / wf_.ecut();
          const double* g2i_ptr = wf_.sd(ispin,ikp)->basis().g2i_ptr();
          double* coeff = (double*) wf_.sd(ispin,ikp)->c().valptr();
          const double* dcoeff =
            (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = wf_.sd(ispin,ikp)->c().mloc();
          const int nloc = wf_.sd(ispin,ikp)->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &coeff[2*mloc*n];
            const double* dc = &dcoeff[2*mloc*n];
            for ( int i = 0; i < mloc; i++ )
            {
              const double g2i = g2i_ptr[i];
              const double dt2bye = ( g2i == 0.0 ? g2i_prec :
                ( g2i < g2i_prec ) ? g2i : g2i_prec );
              c[2*i] -= dt2bye * dc[2*i];
              c[2*i+1] -= dt2bye * dc[2*i+1];
            }
          }
          
          tmap_["gram"].start();
          wf_.sd(ispin,ikp)->gram();
          tmap_["gram"].stop();
        }
      }
    }
  }
}
