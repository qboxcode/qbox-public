////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.C,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#include "SDWavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SDWavefunctionStepper::SDWavefunctionStepper(Sample& s, TimerMap& tmap) : 
  WavefunctionStepper(s,tmap)
{
  dt_ = s_.ctrl.dt;
  const double emass = s_.ctrl.emass;
  dt2bye_ = (emass == 0.0) ? 0.5 / wf_.ecut() : dt_*dt_/emass;
  
  // divide dt2bye by facs coefficient if stress == ON
  if ( s_.ctrl.stress == "ON" )
    dt2bye_ /= s_.ctrl.facs;
}

////////////////////////////////////////////////////////////////////////////////
void SDWavefunctionStepper::update(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      if ( wf_.sd(ispin,ikp) != 0 )
      {
        if ( wf_.sdcontext(ispin,ikp)->active() )
        {
          // c = c - dt2bye * hpsi
          tmap_["update_psi"].start();
          wf_.sd(ispin,ikp)->c().axpy(-dt2bye_,dwf.sd(ispin,ikp)->c());
          tmap_["update_psi"].stop();
          tmap_["gram"].start();
          wf_.sd(ispin,ikp)->gram();
          tmap_["gram"].stop();
        }
      }
    }
  }
}
