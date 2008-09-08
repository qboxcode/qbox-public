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
// $Id: SampleStepper.h,v 1.13 2008-09-08 15:56:19 fgygi Exp $

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
