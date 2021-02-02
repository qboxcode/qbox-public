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
// ANDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ANDIONICSTEPPER_H
#define ANDIONICSTEPPER_H

#include "IonicStepper.h"
#include "AndersonMixer.h"

class ANDIonicStepper : public IonicStepper
{
  private:

  AndersonMixer mixer_;
  double em_;

  public:

  ANDIonicStepper(Sample& s);
  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};
#endif
