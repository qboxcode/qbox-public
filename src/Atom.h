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
// Atom.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ATOM_H
#define ATOM_H

#include "D3vector.h"
#include <string>

class Atom
{
  private:

  std::string name_;
  std::string species_;
  D3vector position_;
  D3vector velocity_;

  public:

  Atom (std::string name, std::string species,
        D3vector position, D3vector velocity);
  std::string name(void) { return name_; };
  std::string species(void) { return species_; };
  D3vector position(void) { return position_; };
  D3vector velocity(void) { return velocity_; };
  void set_position(D3vector p) { position_ = p; };
  void set_velocity(D3vector v) { velocity_ = v; };
  void block(void) { velocity_ = D3vector(0.0,0.0,0.0); };
};

std::ostream& operator << ( std::ostream &os, Atom &a );
#endif
