////////////////////////////////////////////////////////////////////////////////
//
// MDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDWavefunctionStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#ifndef MDWAVEFUNCTIONSTEPPER_H
#define MDWAVEFUNCTIONSTEPPER_H

class Sample;
class Wavefunction;
#include "WavefunctionStepper.h"
#include "Timer.h"
#include <map>
#include <string>
using namespace std;
typedef map<string,Timer> TimerMap;

class MDWavefunctionStepper : public WavefunctionStepper
{
  private:
  Sample& s_;
  Wavefunction& wf_;
  double dt_, dt2bye_;
  TimerMap& tmap_;
  double ekin_ep_, ekin_em_;
  double ekin_eh(void);

  public:

  void update(Wavefunction& dwf);
  void stoermer_start(Wavefunction& dwf);
  void stoermer_end(Wavefunction& dwf);
  double ekin(void) const { return 0.5*(ekin_ep_ + ekin_em_); }

  MDWavefunctionStepper(Sample& s, TimerMap& tmap);
  ~MDWavefunctionStepper() {};
};
#endif
