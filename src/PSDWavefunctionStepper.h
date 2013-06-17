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
// PSDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PSDWAVEFUNCTIONSTEPPER_H
#define PSDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"
class Preconditioner;
class EnergyFunctional;

class PSDWavefunctionStepper : public WavefunctionStepper
{
  private:

  Preconditioner& prec_;

  public:

  void update(Wavefunction& dwf);

  PSDWavefunctionStepper(Wavefunction& wf, Preconditioner& prec,
    TimerMap& tmap);
  ~PSDWavefunctionStepper();
};
#endif
