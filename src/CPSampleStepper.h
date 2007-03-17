////////////////////////////////////////////////////////////////////////////////
//
// CPSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.h,v 1.4 2007-03-17 01:14:00 fgygi Exp $

#ifndef CPSAMPLESTEPPER_H
#define CPSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "ChargeDensity.h"
#include "Sample.h"
#include "Wavefunction.h"
class MDWavefunctionStepper;
class MDIonicStepper;

class CPSampleStepper : public SampleStepper
{
  private:
  
  ChargeDensity cd_;
  EnergyFunctional ef_;
  Wavefunction dwf;
  Wavefunction* wfv;
  
  MDWavefunctionStepper* mdwf_stepper;
  MDIonicStepper* mdionic_stepper;

  // Do not allow construction of CPSampleStepper unrelated to a Sample
  CPSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);

  CPSampleStepper(Sample& s);
  ~CPSampleStepper();
};
#endif
