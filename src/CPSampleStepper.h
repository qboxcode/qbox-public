////////////////////////////////////////////////////////////////////////////////
//
// CPSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.h,v 1.2 2004-03-11 21:52:31 fgygi Exp $

#ifndef CPSAMPLESTEPPER_H
#define CPSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "Sample.h"
#include "Wavefunction.h"
class EnergyFunctional;
class MDWavefunctionStepper;
class MDIonicStepper;
using namespace std;

class CPSampleStepper : public SampleStepper
{
  private:
  
  EnergyFunctional& ef_;
  Wavefunction dwf;
  Wavefunction* wfv;
  
  MDWavefunctionStepper* mdwf_stepper;
  MDIonicStepper* mdionic_stepper;

  // Do not allow construction of CPSampleStepper unrelated to a Sample
  CPSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);

  CPSampleStepper(Sample& s, EnergyFunctional& ef);
  ~CPSampleStepper();
};
#endif
