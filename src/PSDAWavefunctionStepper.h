////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.h,v 1.3 2004-03-11 21:52:31 fgygi Exp $

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
class Preconditioner;

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  Wavefunction wf_last_, dwf_last_;

  // Anderson acceleration flag
  bool extrapolate_;
  
  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) { extrapolate_ = false; }

  PSDAWavefunctionStepper(Sample& s, Preconditioner& p, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif
