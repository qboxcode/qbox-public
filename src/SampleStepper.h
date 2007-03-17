////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.h,v 1.10 2007-03-17 01:14:00 fgygi Exp $

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "Timer.h"
#include <map>
#include <string>
#include <vector>
#include <valarray>

class Sample;
typedef std::map<std::string,Timer> TimerMap;

class SampleStepper
{
  protected:
  
  Sample& s_;
  
  std::vector<std::vector<double> > fion;
  std::valarray<double> sigma_eks, sigma_kin, sigma_ext, sigma;

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
