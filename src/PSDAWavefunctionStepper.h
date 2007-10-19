////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.h,v 1.5 2007-10-19 16:24:04 fgygi Exp $

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
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

  PSDAWavefunctionStepper(Wavefunction& wf, Preconditioner& p, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif
