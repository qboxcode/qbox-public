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
// PairExtForce.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PairExtForce.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void PairExtForce::setup(const AtomSet& atoms)
{
  // find positions in tau array corresponding to atoms name1 and name2
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  is2_ = atoms.is(name2_);
  ia2_ = atoms.ia(name2_);
  assert(is1_>=0);
  assert(ia1_>=0);
  assert(is2_>=0);
  assert(ia2_>=0);
}

////////////////////////////////////////////////////////////////////////////////
double PairExtForce::energy(const vector<vector<double> > &r,
                            vector<vector<double> > &f)
{
  // unit vector of direction joining atom1 and atom2
  const double* pr1 = &r[is1_][3*ia1_];
  const double* pr2 = &r[is2_][3*ia2_];
  D3vector r1(pr1);
  D3vector r2(pr2);

  // A positive value of force_ corresponds to a repulsive force
  // i.e. a positive force tends to increase the distance
  D3vector e12 = normalized(r2-r1);
  f[is1_][3*ia1_+0] -= force_ * e12.x;
  f[is1_][3*ia1_+1] -= force_ * e12.y;
  f[is1_][3*ia1_+2] -= force_ * e12.z;

  f[is2_][3*ia2_+0] += force_ * e12.x;
  f[is2_][3*ia2_+1] += force_ * e12.y;
  f[is2_][3*ia2_+2] += force_ * e12.z;

  return  -2 * force_ * length(r2-r1);
}

////////////////////////////////////////////////////////////////////////////////
ostream& PairExtForce::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <extforce name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"" << name1_ << " " << name2_ << "\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "            value=\"" << setprecision(6) << force_ << "\"/>";
  return os;
}
