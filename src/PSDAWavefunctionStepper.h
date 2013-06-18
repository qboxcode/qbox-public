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
// PSDAWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PSDAWAVEFUNCTIONSTEPPER_H
#define PSDAWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
class Preconditioner;

class PSDAWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;
  Wavefunction wf_last_, dwf_last_;

  // Anderson acceleration flag
  bool extrapolate_;

  public:

  void update(Wavefunction& dwf);
  virtual void preprocess(void) { extrapolate_ = false; }

  PSDAWavefunctionStepper(Wavefunction& wf, Preconditioner& prec,
    TimerMap& tmap);
  ~PSDAWavefunctionStepper();
};
#endif
