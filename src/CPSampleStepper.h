////////////////////////////////////////////////////////////////////////////////
//
// CPSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

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
  
  Wavefunction dwf;
  Wavefunction* wfv;
  vector<vector<double> > fion;
  
  MDWavefunctionStepper* mdwf_stepper;
  MDIonicStepper* mdionic_stepper;

  // Do not allow construction of CPSampleStepper unrelated to a Sample
  CPSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(EnergyFunctional& e, int niter);

  CPSampleStepper(Sample& s);
  //~CPSampleStepper();
};
#endif
