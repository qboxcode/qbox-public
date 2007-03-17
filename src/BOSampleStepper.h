////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.h,v 1.5 2007-03-17 01:14:00 fgygi Exp $

#ifndef BOSAMPLESTEPPER_H
#define BOSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"

class WavefunctionStepper;
class IonicStepper;

class BOSampleStepper : public SampleStepper
{
  private:
  
  Wavefunction dwf;
  Wavefunction* wfv;
  int nitscf_;
  int nite_;
  ChargeDensity cd_;
  EnergyFunctional ef_;
  
  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(int niter);

  BOSampleStepper(Sample& s, int nitscf, int nite);
  ~BOSampleStepper();
};
#endif
