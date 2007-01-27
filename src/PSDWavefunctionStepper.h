////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.4 2007-01-27 23:46:31 fgygi Exp $

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
