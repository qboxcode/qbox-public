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
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BOSAMPLESTEPPER_H
#define BOSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "ChargeDensity.h"
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

  bool initial_atomic_density;
  bool initial_density;
  bool first_step;  // true if step() function has never been called before

  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;

  void step(int niter);
  void initialize_density(void);         // initialize density by atomic density
  void initialize_density(const vector<vector<double> >& rhor); // initialize density by given density

  EnergyFunctional& ef(void) { return ef_; }
  ChargeDensity& cd(void) { return cd_; }

  BOSampleStepper(Sample& s, int nitscf, int nite);
  BOSampleStepper(Sample& s, const vector<vector<double> > * rhor_initial,
                  int nitscf, int nite);
  ~BOSampleStepper();
};
#endif
