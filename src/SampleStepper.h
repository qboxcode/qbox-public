////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.h,v 1.9 2005-09-16 23:08:11 fgygi Exp $

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "Sample.h"
#include "Timer.h"
#include <map>
#include <string>
#include <valarray>
using namespace std;

typedef map<string,Timer> TimerMap;

class SampleStepper
{
  protected:
  
  Sample& s_;
  
  vector<vector<double> > fion;
  valarray<double> sigma_eks, sigma_kin, sigma_ext, sigma;

  // Do not allow construction of SampleStepper unrelated to a Sample
  SampleStepper(void);

  public:

  mutable TimerMap tmap;
  
  virtual void step(int niter) = 0;
  void print_stress(void);
  void compute_sigma(void); // compute kinetic contribution to stress

  SampleStepper(Sample& s);
  virtual ~SampleStepper(void);
};
#endif
