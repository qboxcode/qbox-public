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
// SDAIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SDAIONICSTEPPER_H
#define SDAIONICSTEPPER_H

#include "IonicStepper.h"
#include "LineMinimizer.h"
#include <vector>

class SDAIonicStepper : public IonicStepper
{
  private:

  bool first_step_;
  std::vector<std::vector< double> > rc_;
  std::vector<std::vector< double> > pc_;
  std::vector<std::vector< double> > fc_;
  double ec_, fpc_;
  double alpha_, sigma1_, sigma2_;
  LineMinimizer linmin_;

  public:

  SDAIonicStepper(Sample& s) : IonicStepper(s), first_step_(true),
  sigma1_(0.1), sigma2_(0.3)
  {
    linmin_.set_sigma1(sigma1_);
    linmin_.set_sigma2(sigma2_);
  }

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};

#endif
