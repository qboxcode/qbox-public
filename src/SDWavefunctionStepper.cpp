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
// SDWavefunctionStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SDWavefunctionStepper::SDWavefunctionStepper(Wavefunction& wf, double alpha,
  TimerMap& tmap) :
  alpha_(alpha), WavefunctionStepper(wf,tmap)
{
  tmap_["sd_update_wf"].reset();
  tmap_["gram"].reset();
}

////////////////////////////////////////////////////////////////////////////////
void SDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int isp_loc = 0; isp_loc < wf_.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf_.nkp_loc(); ++ikp_loc )
    {
      // c = c - dt2bye * hpsi
      tmap_["sd_update_wf"].start();
      wf_.sd(isp_loc,ikp_loc)->c().axpy(-alpha_,dwf.sd(isp_loc,ikp_loc)->c());
      tmap_["sd_update_wf"].stop();
      tmap_["gram"].start();
      wf_.sd(isp_loc,ikp_loc)->gram();
      tmap_["gram"].stop();
    }
  }
}
