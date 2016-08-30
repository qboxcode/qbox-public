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

#include <vector>
#include <string>
#include "Sample.h"
#include "ChargeDensity.h"

class Sample;

class ExternalPotential
{
  private:

  Sample& s_;
  int n_[3];                         // real space grid size in 3 dimensions
  double ecut_;
  double amplitude_;
  vector<double> vext_r_;            // vext in real space
  std::string filename_;             // filename for external potential

  public:

  ExternalPotential(Sample& s,std::string name): s_(s),
    filename_(name), ecut_(0.0), amplitude_(1.0) {}
  ~ExternalPotential() {}

  int n(int i) const { return n_[i]; }
  double ecut(void) const { return ecut_; }
  double amplitude(void) const { return amplitude_; }
  std::string filename(void) const { return filename_; }
  double v(size_t i) const { return amplitude_ * vext_r_[i]; }
  void update(const ChargeDensity& cd);
  void set_amplitude(double a) { amplitude_ = a; }
};
#endif
