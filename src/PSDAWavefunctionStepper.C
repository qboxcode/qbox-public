////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.C,v 1.4 2004-02-04 19:55:16 fgygi Exp $

#include "PSDAWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDAWavefunctionStepper::PSDAWavefunctionStepper(Sample& s, TimerMap& tmap) : 
  WavefunctionStepper(s,tmap), wf_last_(s.wf), dwf_last_(s.wf), 
  extrapolate_(false)
{
  dt_ = s_.ctrl.dt;
  const double emass = s_.ctrl.emass;
  dt2bye_ = (emass == 0.0) ? 0.5 / wf_.ecut() : dt_*dt_/emass;
  
  // divide dt2bye by facs coefficient if stress == ON
  if ( s_.ctrl.stress == "ON" )
    dt2bye_ /= s_.ctrl.facs;
}

////////////////////////////////////////////////////////////////////////////////
void PSDAWavefunctionStepper::update(Wavefunction& dwf)
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
 
            // dwf.sd->c() now contains the descent direction (HV-VA)
 
            // Apply preconditioner K and store dt2bye*K(HV-VA) in dwf
 
            const double g2i_prec = s_.ctrl.ecutprec > 0.0 ?
                                    0.5 / s_.ctrl.ecutprec :
                                    0.5 / wf_.ecut();
            const double* g2i_ptr = wf_.sd(ispin,ikp)->basis().g2i_ptr();
            double* c = (double*) wf_.sd(ispin,ikp)->c().valptr();
            double* c_last = (double*) wf_last_.sd(ispin,ikp)->c().valptr();
            double* dc = (double*) dwf.sd(ispin,ikp)->c().valptr();
            double* dc_last = (double*) dwf_last_.sd(ispin,ikp)->c().valptr();
            const int mloc = wf_.sd(ispin,ikp)->c().mloc();
            const int nloc = wf_.sd(ispin,ikp)->c().nloc();
            
            // next line: add enhancement factor to descent direction
            // since there is no instability of the Anderson iteration
            // This improves convergence in most cases
            const double psda_enhancement_factor = 2.0;
            for ( int n = 0; n < nloc; n++ )
            {
              // note: double mloc length for complex<double> indices
              double* dcn = &dc[2*mloc*n];
              for ( int i = 0; i < mloc; i++ )
              {
                const double g2i = g2i_ptr[i];
                const double dt2bye = ( g2i == 0.0 ? g2i_prec :
                  ( g2i < g2i_prec ) ? g2i : g2i_prec );
                const double f0 = -psda_enhancement_factor*dt2bye * dcn[2*i];
                const double f1 = -psda_enhancement_factor*dt2bye * dcn[2*i+1];
                dcn[2*i] = f0;
                dcn[2*i+1] = f1;
              }
            }
            // dwf now contains the preconditioned descent
            // direction -dt2bye*K(HV-VA)
 
            // Anderson extrapolation
            if ( extrapolate_ )
            {
              double theta = 0.0;
              double a = 0.0, b = 0.0;
              for ( int i = 0; i < 2*mloc*nloc; i++ )
              {
                const double f = dc[i];
                const double delta_f = f - dc_last[i];
 
                // accumulate partial sums of a and b
                // a = delta_F * F
 
                a += f * delta_f;
                b += delta_f * delta_f;
              }
              // correct for double counting of asum and bsum on first row
              // factor 2.0: G and -G
              a *= 2.0;
              b *= 2.0;
              if ( wf_.sdcontext(ispin,ikp)->myrow() == 0 )
              {
                for ( int n = 0; n < nloc; n++ )
                {
                  const int i = 2*mloc*n;
                  const double f0 = dc[i];
                  const double f1 = dc[i+1];
                  const double delta_f0 = f0 - dc_last[i];
                  const double delta_f1 = f1 - dc_last[i+1];
                  a -= f0 * delta_f0 + f1 * delta_f1;
                  b -= delta_f0 * delta_f0 + delta_f1 * delta_f1;
                }
              }
 
              // a and b contain the partial sums of a and b
              double tmpvec[2] = { a, b };
              wf_.sdcontext(ispin,ikp)->dsum(2,1,&tmpvec[0],1);
              a = tmpvec[0];
              b = tmpvec[1];
 
              // compute theta = - a / b
              if ( b != 0.0 )
                theta = - a / b;
 
              if ( wf_.sdcontext(ispin,ikp)->onpe0() )
                cout << "  <!-- Anderson extrapolation: theta=" << theta;
 
              if ( theta < -2.0 )
              {
                // cancel extrapolation for the current
                // update by setting theta = 0.0
                theta = 0.0;
              }
 
              theta = max(-0.5,theta);
              theta = min(3.0,theta);
 
              if ( wf_.sdcontext(ispin,ikp)->onpe0() )
                cout <<" (" << theta << ")" << " -->"<< endl;
                
              // extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ )
              {
                // x_bar = x_ + theta * ( x_ - xlast_ ) (store in x_)
                const double x = c[i];
                const double xlast = c_last[i];
                const double xbar = x + theta * ( x - xlast );
 
                // f_bar = f + theta * ( f - flast ) (store in f)
                const double f = dc[i];
                const double flast = dc_last[i];
                const double fbar = f + theta * ( f - flast );
 
                c[i] = xbar + fbar;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            else
            {
              // no extrapolation
              for ( int i = 0; i < 2*mloc*nloc; i++ )
              {
                // x_ = x_ + f_
                const double x = c[i];
                const double f = dc[i];
 
                c[i] = x + f;
                c_last[i] = x;
                dc_last[i] = f;
              }
            }
            extrapolate_ = true;
 
            tmap_["riccati"].start();
            wf_.sd(ispin,ikp)->riccati(*wf_last_.sd(ispin,ikp));
            tmap_["riccati"].stop();
          }
          else
          {
            // not implemented in the complex case
            assert(false);
          }
        }
      }
    }
  }
}
