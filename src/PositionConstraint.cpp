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
// PositionConstraint.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PositionConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void PositionConstraint::setup(const AtomSet& atoms)
{
  // find position in tau array corresponding to atom name1
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  assert(is1_>=0);
  assert(ia1_>=0);
}

////////////////////////////////////////////////////////////////////////////////
void PositionConstraint::update(double dt)
{
  // nothing to update
}

////////////////////////////////////////////////////////////////////////////////
bool PositionConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  double* pr1p  = &rp[is1_][3*ia1_];
  D3vector r1p(pr1p);

  double sigma = length(r1p-r1);

  if ( sigma < tol_ ) return true;

  pr1p[0] = pr1[0];
  pr1p[1] = pr1[1];
  pr1p[2] = pr1[2];

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool PositionConstraint::enforce_v(const vector<vector<double> > &r0,
vector<vector<double> > &v0) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  double* pv1 = &v0[is1_][3*ia1_];
  D3vector v1(pv1);

  const double err = length(v1);
  if ( err < tol_ ) return true;

  pv1[0] = 0.0;
  pv1[1] = 0.0;
  pv1[2] = 0.0;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void PositionConstraint::compute_force(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  const double* pf1 = &f[is1_][3*ia1_];
  D3vector f1(pf1);

  force_ = length(f1);
}

////////////////////////////////////////////////////////////////////////////////
ostream& PositionConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type() << "\">";
  os << " <atoms> " << name1_ << " </atoms>\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  <velocity> " << setprecision(8) << 0 << " </velocity>";
  os << " <weight> " << setprecision(8) << weight_ << " </weight>\n";
  os << "  <value> " << setprecision(8) << 0 << " </value>";
  os << " <force> " << setprecision(8) << force_ << " </force>\n";
  os << " </constraint>";
  return os;
}
