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
// PSDAWavefunctionStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PSDAWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDAWavefunctionStepper::PSDAWavefunctionStepper(Wavefunction& wf,
  Preconditioner& prec, TimerMap& tmap) : prec_(prec),
  WavefunctionStepper(wf,tmap), wf_last_(wf), dwf_last_(wf),
  extrapolate_(false)
{
  tmap_["psda_residual"].reset();
  tmap_["psda_prec"].reset();
  tmap_["psda_update_wf"].reset();
  tmap_["gram"].reset();
  tmap_["lowdin"].reset();
  tmap_["ortho_align"].reset();
  tmap_["riccati"].reset();
}

////////////////////////////////////////////////////////////////////////////////
PSDAWavefunctionStepper::~PSDAWavefunctionStepper(void)
{}

////////////////////////////////////////////////////////////////////////////////
void PSDAWavefunctionStepper::update(Wavefunction& dwf)
{
  tmap_["psda_residual"].start();
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      // compute A = V^T H V  and descent direction HV - VA
      SlaterDet* sd = wf_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd = dwf.sd(isp_loc,ikp_loc);
      if ( sd->basis().real() )
      {
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
      else
      {
        ComplexMatrix& c = sd->c();
        ComplexMatrix& cp = dsd->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        a.gemm('c','n',1.0,c,cp,0.0);
        // cp = cp - c * a
        cp.gemm('n','n',-1.0,c,a,1.0);
      }
    }
  }
  tmap_["psda_residual"].stop();

  // dwf.sd->c() now contains the descent direction (HV-VA) (residual)
  // update the preconditioner
  prec_.update(wf_);

  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      tmap_["psda_prec"].start();
      SlaterDet* sd = wf_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd = dwf.sd(isp_loc,ikp_loc);
      SlaterDet* sd_last = wf_last_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd_last = dwf_last_.sd(isp_loc,ikp_loc);
      // Apply preconditioner K and store -K(HV-VA) in dwf
      double* c = (double*) sd->c().valptr();
      double* c_last = (double*) sd_last->c().valptr();
      double* dc = (double*) dsd->c().valptr();
      double* dc_last = (double*) dsd_last->c().valptr();
      const int mloc = sd->c().mloc();
      const int ngwl = sd->basis().localsize();
      const int nloc = sd->c().nloc();

      for ( int n = 0; n < nloc; n++ )
      {
        // note: double mloc length for complex<double> indices
        double* dcn = &dc[2*mloc*n];
        for ( int i = 0; i < ngwl; i++ )
        {
          const double fac = prec_.diag(isp_loc,ikp_loc,n,i);
          const double f0 = -fac * dcn[2*i];
          const double f1 = -fac * dcn[2*i+1];
          dcn[2*i] = f0;
          dcn[2*i+1] = f1;
        }
      }
      tmap_["psda_prec"].stop();

      // dwf now contains the preconditioned descent
      // direction -K(HV-VA)

      tmap_["psda_update_wf"].start();
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

        if ( sd->basis().real() )
        {
          // correct for double counting of asum and bsum on first row
          // factor 2.0: G and -G
          a *= 2.0;
          b *= 2.0;
          if ( wf_.sd_context().myrow() == 0 )
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
        }

        // a and b contain the partial sums of a and b
        double tmpvec[2] = { a, b };
        wf_.sd_context().dsum(2,1,&tmpvec[0],2);
        a = tmpvec[0];
        b = tmpvec[1];

        // compute theta = - a / b
        if ( b != 0.0 )
          theta = - a / b;

        if ( theta < -1.0 )
        {
          theta = 0.0;
        }

        theta = min(2.0,theta);

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
      tmap_["psda_update_wf"].stop();

      enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };
      //const ortho_type ortho = GRAM;
      //const ortho_type ortho = LOWDIN;
      const ortho_type ortho = ORTHO_ALIGN;

      switch ( ortho )
      {
        case GRAM:
          tmap_["gram"].start();
          sd->gram();
          tmap_["gram"].stop();
          break;

        case LOWDIN:
          tmap_["lowdin"].start();
          sd->lowdin();
          tmap_["lowdin"].stop();
          break;

        case ORTHO_ALIGN:
          tmap_["ortho_align"].start();
          sd->ortho_align(*sd_last);
          tmap_["ortho_align"].stop();
          break;

        case RICCATI:
          tmap_["riccati"].start();
          sd->riccati(*sd_last);
          tmap_["riccati"].stop();
          break;
      }
    } // ikp_loc
  } // isp_loc
  extrapolate_ = true;
}
