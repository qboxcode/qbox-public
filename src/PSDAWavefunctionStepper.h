////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.h,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:

  Wavefunction wf_last_, dwf_last_;
  double dt_, dt2bye_;

  // Anderson acceleration flag
  bool extrapolate_;
  
  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) { extrapolate_ = false; }

  PSDAWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif
