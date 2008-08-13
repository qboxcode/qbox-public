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
// WavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.h,v 1.9 2008-08-13 06:39:43 fgygi Exp $

#ifndef WAVEFUNCTIONSTEPPER_H
#define WAVEFUNCTIONSTEPPER_H
#include "Timer.h"
#include <map>
#include <string>

typedef std::map<std::string,Timer> TimerMap;
class Wavefunction;

class WavefunctionStepper
{
  private:

  protected:
  Wavefunction& wf_;
  TimerMap& tmap_;

  public:

  virtual void update(Wavefunction& dwf) = 0;
  virtual void preprocess(void) {}
  virtual void postprocess(void) {}

  WavefunctionStepper(Wavefunction& wf, TimerMap& tmap) :
  wf_(wf), tmap_(tmap)
  {}
  virtual ~WavefunctionStepper() {}
};
#endif
