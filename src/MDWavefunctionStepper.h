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
// MDWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MDWAVEFUNCTIONSTEPPER_H
#define MDWAVEFUNCTIONSTEPPER_H

#include "WavefunctionStepper.h"

class MDWavefunctionStepper : public WavefunctionStepper
{
  private:

  double dt_;
  double dt2bye_;
  Wavefunction *wfv_;

  double ekin_ep_, ekin_em_;
  double ekin_eh(void);

  public:

  void update(Wavefunction& dwf);
  void compute_wfm(Wavefunction& dwf);
  void compute_wfv(Wavefunction& dwf);
  double ekin(void) const { return 0.5*(ekin_ep_ + ekin_em_); }

  MDWavefunctionStepper(Wavefunction& wf, Wavefunction* wfv,
    double dt, double dt2bye, TimerMap& tmap);
  ~MDWavefunctionStepper() {};
};
#endif
