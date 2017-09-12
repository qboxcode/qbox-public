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
// ExtForce.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EXTFORCE_H
#define EXTFORCE_H

#include <string>
#include <vector>
#include <cassert>
#include "D3vector.h"

class AtomSet;

class ExtForce
{
  protected:

  std::string name_;               // extforce name
  std::vector<std::string> names_; // names of atoms involved in the extforce

  public:

  virtual ~ExtForce(void){}
  virtual std::string type(void) const = 0;
  virtual void set_value(std::vector<double> f) = 0;
  virtual double energy(const std::vector<std::vector<double> > &r,
                        std::vector<std::vector<double> > &f) = 0;
  virtual void setup(const AtomSet& atoms) = 0;
  virtual std::ostream& print(std::ostream &os) = 0;
  virtual D3vector sum_contrib(void) = 0;
  std::string name(void) const { return name_; }
  virtual double magnitude(void) = 0;
  std::string names(int i) const
  {
    assert( i >= 0 && i < names_.size() );
    return names_[i];
  }
};
std::ostream& operator << ( std::ostream &os, ExtForce &xf );
#endif
