////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.h,v 1.2 2004-03-11 21:52:32 fgygi Exp $

#ifndef BOSAMPLESTEPPER_H
#define BOSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "Sample.h"
#include "Wavefunction.h"
class EnergyFunctional;
class WavefunctionStepper;
class IonicStepper;
using namespace std;

class BOSampleStepper : public SampleStepper
{
  private:
  
  Wavefunction dwf;
  Wavefunction* wfv;
  int nite_;
  EnergyFunctional& ef_;
  
  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);

  BOSampleStepper(Sample& s, EnergyFunctional& ef, int nite);
  //~BOSampleStepper();
};
#endif
