////////////////////////////////////////////////////////////////////////////////
//
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDAWavefunctionStepper.h,v 1.1 2003-11-27 01:30:31 fgygi Exp $

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

class Sample;
#include "Wavefunction.h"
#include "WavefunctionStepper.h"
#include "Timer.h"
#include <map>
#include <string>
#include <valarray>
using namespace std;
typedef map<string,Timer> TimerMap;

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:
  Sample& s_;
  Wavefunction& wf_;
  Wavefunction wf_last_, dwf_last_;
  double dt_, dt2bye_;
  TimerMap& tmap_;

  // Anderson acceleration flag
  bool extrapolate_;
  
  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) { extrapolate_ = false; }

  PSDAWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~PSDAWavefunctionStepper() {};
};
#endif
