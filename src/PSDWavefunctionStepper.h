////////////////////////////////////////////////////////////////////////////////
//
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: PSDWavefunctionStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#ifndef PSDWAVEFUNCTIONSTEPPER_H
#define PSDWAVEFUNCTIONSTEPPER_H

class Sample;
class Wavefunction;
#include "WavefunctionStepper.h"
#include "Timer.h"
#include <map>
#include <string>
using namespace std;
typedef map<string,Timer> TimerMap;

class PSDWavefunctionStepper : public WavefunctionStepper
{
  private:
  Sample& s_;
  Wavefunction& wf_;
  double dt_, dt2bye_;
  TimerMap& tmap_;

  public:

  void update(Wavefunction& dwf);

  PSDWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~PSDWavefunctionStepper() {};
};
#endif
