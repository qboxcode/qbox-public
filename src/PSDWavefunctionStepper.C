////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.C,v 1.6 2004-11-10 22:35:23 fgygi Exp $

#include "PSDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDWavefunctionStepper::PSDWavefunctionStepper(Sample& s, 
  Preconditioner& p, TimerMap& tmap) : 
  WavefunctionStepper(s,tmap), prec_(p)
{}

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
 
          tmap_["psd_residual"].start();
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
          tmap_["psd_residual"].stop();
 
          // dwf.sd->c() now contains the descent direction (HV-VA)
 
          tmap_["psd_update_wf"].start();
          const valarray<double>& diag = prec_.diag(ispin,ikp);
          
          double* coeff = (double*) wf_.sd(ispin,ikp)->c().valptr();
          const double* dcoeff =
            (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
          const int mloc = wf_.sd(ispin,ikp)->c().mloc();
          const int ngwl = wf_.sd(ispin,ikp)->basis().localsize();
          const int nloc = wf_.sd(ispin,ikp)->c().nloc();
          for ( int n = 0; n < nloc; n++ )
          {
            // note: double mloc length for complex<double> indices
            double* c = &coeff[2*mloc*n];
            const double* dc = &dcoeff[2*mloc*n];
            // loop to ngwl only since diag[i] is not defined on [0:mloc-1]
            for ( int i = 0; i < ngwl; i++ )
            {
              const double fac = diag[i];
              const double delta_re = fac * dc[2*i];
              const double delta_im = fac * dc[2*i+1];
              c[2*i]   -= delta_re;
              c[2*i+1] -= delta_im;
            }
          }
          tmap_["psd_update_wf"].stop();
          
          tmap_["gram"].start();
          wf_.sd(ispin,ikp)->gram();
          tmap_["gram"].stop();
        }
      }
    }
  }
}
