////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

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
  vector<vector<double> > tau0,taum,vel,fion;
  int nite_;
  
  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(EnergyFunctional& e, int niter);

  BOSampleStepper(Sample& s, int nite);
  //~BOSampleStepper();
};
#endif
