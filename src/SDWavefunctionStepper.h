////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.h,v 1.3 2007-01-27 23:46:31 fgygi Exp $

#ifndef SDWAVEFUNCTIONSTEPPER_H
#define SDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class SDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double alpha_;

  public:

  void update(Wavefunction& dwf);

  SDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap);
  ~SDWavefunctionStepper() {};
};
#endif
