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
// SDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SDWAVEFUNCTIONSTEPPER_H
#define SDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class SDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double alpha_;

  public:

  void update(Wavefunction& dwf);

  SDWavefunctionStepper(Wavefunction& wf, double alpha, TimerMap& tmap);
  ~SDWavefunctionStepper() {};
};
#endif
