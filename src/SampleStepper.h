////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.h,v 1.7 2003-11-21 20:01:47 fgygi Exp $

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "Sample.h"
#include "Timer.h"
class EnergyFunctional;
#include <map>
#include <string>
using namespace std;

typedef map<string,Timer> TimerMap;

class SampleStepper
{
  protected:
  
  Sample& s_;

  // Do not allow construction of SampleStepper unrelated to a Sample
  SampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  virtual void step(EnergyFunctional& e, int niter) = 0;

  SampleStepper(Sample& s);
  virtual ~SampleStepper(void);
};
#endif
