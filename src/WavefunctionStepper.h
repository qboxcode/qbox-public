////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.5 2004-03-11 21:52:31 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
#include "Sample.h"
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
  Sample& s_;
  Wavefunction& wf_;
  TimerMap& tmap_;
  
  public:

  virtual void update(Wavefunction& dwf) = 0;
  virtual void preprocess(void) {}
  virtual void postprocess(void) {}

  WavefunctionStepper(Sample& s, TimerMap& tmap) : 
  s_(s), wf_(s.wf), tmap_(tmap)
  {}
  virtual ~WavefunctionStepper() {}
};
#endif
