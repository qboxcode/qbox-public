////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#ifndef PSDWAVEFUNCTIONSTEPPER_H
#define PSDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class PSDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double dt_, dt2bye_;

  public:

  void update(Wavefunction& dwf);

  PSDWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~PSDWavefunctionStepper() {};
};
#endif
