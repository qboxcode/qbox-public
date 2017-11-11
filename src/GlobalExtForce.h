////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
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
// GlobalExtForce.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef GLOBALEXTFORCE_H
#define GLOBALEXTFORCE_H

#include "ExtForce.h"
#include "D3vector.h"
#include <cassert>

class AtomSet;

class GlobalExtForce: public ExtForce
{
  D3vector force_;
  int na_;

  public:

  GlobalExtForce(std::string name, D3vector force) : force_(force)
  {
    name_ = name;
  }
  ~GlobalExtForce(void) {}

  std::string type(void) const { return "global"; }
  void set_value(std::vector<double> f)
  {
    force_.x = f[0];
    force_.y = f[1];
    force_.z = f[2];
  }

  void setup(const AtomSet& atoms);
  double energy(const std::vector<std::vector<double> > &r,
                std::vector<std::vector<double> > &f);
  D3vector sum_contrib(void) { return force_; }
  double magnitude(void) { return length(force_); }
  std::ostream& print( std::ostream& os );

};
#endif
