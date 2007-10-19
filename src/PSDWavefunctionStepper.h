////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.5 2007-10-19 16:24:04 fgygi Exp $

#ifndef PSDWAVEFUNCTIONSTEPPER_H
#define PSDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
class Preconditioner;

class PSDWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;

  public:

  void update(Wavefunction& dwf);

  PSDWavefunctionStepper(Wavefunction& wf, Preconditioner& p, TimerMap& tmap);
  ~PSDWavefunctionStepper() {};
};
#endif
