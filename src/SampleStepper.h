////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.h,v 1.6 2003-06-11 22:10:11 fgygi Exp $

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include <iostream>
#include <map>
#include <string>
using namespace std;

typedef map<string,Timer> TimerMap;

class SampleStepper
{
  private:
  
  Sample& s_;
  Wavefunction dwf;
  Wavefunction* wfv;
  vector<vector<double> > tau0,taum,vel,fion;
  vector<double> pmass;
  UnitCell dcell;

  // Do not allow construction of SampleStepper unrelated to a Sample
  SampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  void step(EnergyFunctional& e, int niter);

  SampleStepper(Sample& s);
  ~SampleStepper();
};
#endif
