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
// PSDWavefunctionStepper.C
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
{}

////////////////////////////////////////////////////////////////////////////////
PSDWavefunctionStepper::~PSDWavefunctionStepper(void)
{}

////////////////////////////////////////////////////////////////////////////////
void PSDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
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
        ComplexMatrix& c = wf_.sd(ispin,ikp)->c();
        ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
        ComplexMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
        a.gemm('c','n',1.0,c,cp,0.0);
        // cp = cp - c * a
        cp.gemm('n','n',-1.0,c,a,1.0);
      }
      tmap_["psd_residual"].stop();
    }
  }

  // dwf.sd->c() now contains the descent direction (HV-VA)

  // update preconditioner using the residual
  prec_.update(dwf);

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      tmap_["psd_update_wf"].start();
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

        for ( int i = 0; i < ngwl; i++ )
        {
          const double fac = prec_.diag(ispin,ikp,n,i);
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
