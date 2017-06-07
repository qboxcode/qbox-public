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
// PairExtForce.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PAIREXTFORCE_H
#define PAIREXTFORCE_H

#include "ExtForce.h"
#include "D3vector.h"
#include <cassert>

class AtomSet;

class PairExtForce: public ExtForce
{
  std::string name1_, name2_; // names of atoms
  int    ia1_, is1_, ia2_, is2_;
  double force_;

  public:

  PairExtForce(std::string name, std::string name1, std::string name2,
    double force) :
  name1_(name1), name2_(name2), force_(force)
  {
    name_ = name;
    names_.resize(2);
    names_[0] = name1_;
    names_[1] = name2_;
  }
  ~PairExtForce(void) {}

  std::string type(void) const { return "pair"; }
  void set_value(std::vector<double> f)
  {
    assert(f.size()==1);
    force_= f[0];
  }

  void setup(const AtomSet& atoms);
  double energy(const std::vector<std::vector<double> > &r,
                std::vector<std::vector<double> > &f);
  D3vector sum_contrib(void) { return D3vector(0.,0.,0.); }
  double magnitude(void) { return fabs(force_); }
  std::ostream& print( std::ostream& os );

};
#endif
