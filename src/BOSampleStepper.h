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
  int nitscf_;
  int nite_;
  ChargeDensity cd_;
  EnergyFunctional ef_;

  WavefunctionStepper* wf_stepper;
  IonicStepper* ionic_stepper;

  bool update_density_first_;
  bool update_vh_, update_vxc_;

  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;

  void step(int niter);
  // initialize density with sum of atomic densities
  void initialize_density(void);
  void set_update_vh(bool update_vh) { update_vh_ = update_vh; }
  void set_update_vxc(bool update_vxc) { update_vxc_ = update_vxc; }
  void set_update_density_first(bool update_density_first)
    { update_density_first_ = update_density_first; }

  EnergyFunctional& ef(void) { return ef_; }
  ChargeDensity& cd(void) { return cd_; }

  BOSampleStepper(Sample& s, int nitscf, int nite);
  ~BOSampleStepper();
};
#endif
