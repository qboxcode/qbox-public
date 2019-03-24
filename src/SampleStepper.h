////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLESTEPPER_H
#define SAMPLESTEPPER_H

#include "Timer.h"
#include "EnergyFunctional.h"
#include "ChargeDensity.h"
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

  std::string iter_cmd_;
  int iter_cmd_period_;

  // Do not allow construction of SampleStepper unrelated to a Sample
  SampleStepper(void);

  public:

  mutable TimerMap tmap;

  virtual void step(int niter) = 0;
  void print_stress(void);
  void compute_sigma(void); // compute kinetic contribution to stress
  virtual void initialize_density() {}
  void set_iter_cmd(std::string s) { iter_cmd_ = s; }
  void set_iter_cmd_period(int i) { iter_cmd_period_ = i; }

  virtual EnergyFunctional& ef(void) = 0;
  virtual ChargeDensity& cd(void) = 0;

  SampleStepper(Sample& s);
  virtual ~SampleStepper(void);
};
#endif
