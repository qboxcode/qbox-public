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
// PSDWavefunctionStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PSDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Preconditioner.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
PSDWavefunctionStepper::PSDWavefunctionStepper(Wavefunction& wf,
  Preconditioner& prec, TimerMap& tmap) : WavefunctionStepper(wf,tmap),
  prec_(prec)
{
  tmap_["psd_residual"].reset();
  tmap_["psd_update_wf"].reset();
  tmap_["gram"].reset();
}

////////////////////////////////////////////////////////////////////////////////
PSDWavefunctionStepper::~PSDWavefunctionStepper(void)
{}

////////////////////////////////////////////////////////////////////////////////
void PSDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      // compute A = V^T H V  and descent direction HV - VA
      SlaterDet* sd = wf_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd = dwf.sd(isp_loc,ikp_loc);
      tmap_["psd_residual"].start();
      if ( sd->basis().real() )
      {
        // proxy real matrices c, cp
        DoubleMatrix c(sd->c());
        DoubleMatrix cp(dsd->c());

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
        ComplexMatrix& c = sd->c();
        ComplexMatrix& cp = dsd->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        a.gemm('c','n',1.0,c,cp,0.0);
        // cp = cp - c * a
        cp.gemm('n','n',-1.0,c,a,1.0);
      }
      tmap_["psd_residual"].stop();
    }
  }

  // dwf.sd->c() now contains the descent direction (HV-VA)

  // update preconditioner
  prec_.update(wf_);

  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      SlaterDet* sd = wf_.sd(isp_loc,ikp_loc);
      SlaterDet* dsd = dwf.sd(isp_loc,ikp_loc);
      tmap_["psd_update_wf"].start();
      double* coeff = (double*) sd->c().valptr();
      const double* dcoeff = (const double*) dsd->c().cvalptr();
      const int mloc = sd->c().mloc();
      const int ngwl = sd->basis().localsize();
      const int nloc = sd->c().nloc();
      for ( int n = 0; n < nloc; n++ )
      {
        // note: double mloc length for complex<double> indices
        double* c = &coeff[2*mloc*n];
        const double* dc = &dcoeff[2*mloc*n];

        for ( int i = 0; i < ngwl; i++ )
        {
          const double fac = prec_.diag(isp_loc,ikp_loc,n,i);
          const double delta_re = fac * dc[2*i];
          const double delta_im = fac * dc[2*i+1];
          c[2*i]   -= delta_re;
          c[2*i+1] -= delta_im;
        }
      }
      tmap_["psd_update_wf"].stop();

      tmap_["gram"].start();
      sd->gram();
      tmap_["gram"].stop();
    }
  }
}
