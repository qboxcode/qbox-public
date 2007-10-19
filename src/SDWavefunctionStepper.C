////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.C,v 1.6 2007-10-19 17:37:06 fgygi Exp $

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
