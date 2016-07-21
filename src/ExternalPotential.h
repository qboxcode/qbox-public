////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2016 The Regents of the University of California
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
// ExternalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
#ifndef EXTERNALPOTENTIAL_H
#define EXTERNALPOTENTIAL_H

#include "Sample.h"
#include "ChargeDensity.h"
#include <vector>
#include <string>

class ExternalPotential
{
  private:

  Sample& s_;
  int n_[3];
  std::vector<double> vext_r_;
  std::vector<complex<double> > vext_g_;

  public:

  ExternalPotential(Sample& s): s_(s) {}
  ~ExternalPotential() {}
  double v(size_t i) const { return vext_r_[i]; }
  void update(const ChargeDensity& cd);
};
#endif
