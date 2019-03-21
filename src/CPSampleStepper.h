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
// CPSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CPSAMPLESTEPPER_H
#define CPSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "ChargeDensity.h"
#include "Sample.h"
#include "Wavefunction.h"
class MDWavefunctionStepper;
class MDIonicStepper;

class CPSampleStepper : public SampleStepper
{
  private:

  ChargeDensity cd_;
  EnergyFunctional ef_;
  Wavefunction dwf;
  Wavefunction* wfv;

  MDWavefunctionStepper* mdwf_stepper;
  MDIonicStepper* mdionic_stepper;

  // Do not allow construction of CPSampleStepper unrelated to a Sample
  CPSampleStepper(void);

  public:

  mutable TimerMap tmap;

  void step(int niter);
  ChargeDensity& cd(void) { return cd_; }
  EnergyFunctional& ef(void) { return ef_; }

  CPSampleStepper(Sample& s);
  ~CPSampleStepper();
};
#endif
