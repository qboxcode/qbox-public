////////////////////////////////////////////////////////////////////////////////
//
// SDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDWavefunctionStepper.h,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#ifndef SDWAVEFUNCTIONSTEPPER_H
#define SDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class SDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double dt_, dt2bye_;

  public:

  void update(Wavefunction& dwf);

  SDWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~SDWavefunctionStepper() {};
};
#endif
