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
// AtomicExtForce.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "AtomicExtForce.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void AtomicExtForce::setup(const AtomSet& atoms)
{
  // find position in tau array corresponding to atom name1
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  assert(is1_>=0);
  assert(ia1_>=0);
}

////////////////////////////////////////////////////////////////////////////////
double AtomicExtForce::energy(const vector<vector<double> > &r,
                              vector<vector<double> > &f)
{
  f[is1_][3*ia1_+0] += force_.x;
  f[is1_][3*ia1_+1] += force_.y;
  f[is1_][3*ia1_+2] += force_.z;
  return -force_.x*r[is1_][3*ia1_+0] +
          force_.y*r[is1_][3*ia1_+1] +
          force_.z*r[is1_][3*ia1_+2];
}

////////////////////////////////////////////////////////////////////////////////
ostream& AtomicExtForce::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <extforce name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << name1_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "            value=\"" << setprecision(6) << force_.x << " "
     << setprecision(6) << force_.y << " "
     << setprecision(6) << force_.z << "\"/>";
  return os;
}
