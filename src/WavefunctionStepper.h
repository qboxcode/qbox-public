////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.6 2007-01-27 23:46:31 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
#include "Timer.h"
#include <map>
#include <string>
using namespace std;

typedef map<string,Timer> TimerMap;
class Wavefunction;

class WavefunctionStepper
{
  private:
  
  protected:
  Wavefunction& wf_;
  TimerMap& tmap_;
  
  public:

  virtual void update(Wavefunction& dwf) = 0;
  virtual void preprocess(void) {}
  virtual void postprocess(void) {}

  WavefunctionStepper(Wavefunction& wf, TimerMap& tmap) : 
  wf_(wf), tmap_(tmap)
  {}
  virtual ~WavefunctionStepper() {}
};
#endif
