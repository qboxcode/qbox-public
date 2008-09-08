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
// SDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.C,v 1.8 2008-09-08 15:56:19 fgygi Exp $

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
{}

////////////////////////////////////////////////////////////////////////////////
void SDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      // c = c - dt2bye * hpsi
      tmap_["sd_update_wf"].start();
      wf_.sd(ispin,ikp)->c().axpy(-alpha_,dwf.sd(ispin,ikp)->c());
      tmap_["sd_update_wf"].stop();
      tmap_["gram"].start();
      wf_.sd(ispin,ikp)->gram();
      tmap_["gram"].stop();
    }
  }
}
