////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.3 2004-03-11 21:52:32 fgygi Exp $

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
  
  PSDWavefunctionStepper(Sample& s, Preconditioner& p, TimerMap& tmap);
  ~PSDWavefunctionStepper() {};
};
#endif
