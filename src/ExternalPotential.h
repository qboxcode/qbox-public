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
  int n_[3];               // real space grid size in 3 dimensions
                           // read from cube file in cube file mode,
                           // otherwise must be given in constructor
  double ecut_;
  double magnitude_;       // the magnitude of external potential, defined as
                           // the average of its largest 0.1% (absolute) values
  double amplitude_;       // overall scaling factor of external potential
  vector<double> vext_r_;  // vext in real space
  std::string filename_;   // file name for external potential
  std::string fmt_;        // file format: "cube" or "xml"

  public:

  ExternalPotential(Sample& s, std::string name, std::string fmt="xml"):
    s_(s), filename_(name), ecut_(0.0), amplitude_(1.0), magnitude_(0.0),
    fmt_(fmt)
  {
    assert( fmt_ == "cube" || fmt_ == "xml" );
  }
  ~ExternalPotential() {}

  int n(int i) const { return n_[i]; }
  double ecut(void) const { return ecut_; }
  double magnitude(void) const { return magnitude_; }
  double amplitude(void) const { return amplitude_; }
  std::string filename(void) const { return filename_; }
  double v(size_t i) const { return amplitude_ * vext_r_[i]; }
  void update(const ChargeDensity& cd);
  void set_amplitude(double a) { amplitude_ = a; }
  void reverse(void) {amplitude_ *= -1; }
  double compute_eext(const ChargeDensity& cd);
};
#endif
