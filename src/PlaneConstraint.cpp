////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2024 The Regents of the University of California
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
// PlaneConstraint.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PlaneConstraint.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void PlaneConstraint::setup(const AtomSet& atoms)
{
  // find position in tau array corresponding to atom name1
  is1_ = atoms.is(name1_);
  ia1_ = atoms.ia(name1_);
  assert(is1_>=0);
  assert(ia1_>=0);
}

////////////////////////////////////////////////////////////////////////////////
void PlaneConstraint::update(double dt)
{
  distance_ += velocity_ * dt;
}

////////////////////////////////////////////////////////////////////////////////
bool PlaneConstraint::enforce_r(const vector<vector<double> > &r0,
vector<vector<double> > &rp) const
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  double* pr1p  = &rp[is1_][3*ia1_];
  D3vector r1p(pr1p);

  // project r1 onto plane defined by e_ : ex * x + ey * y + ez * z = d
  // pr1p = pr1 - e_ * pr1

  // sigma = distance to plane: sigma = r1p * e_ - distance_
  double sigma = r1p * e_ - distance_;

  if ( fabs(sigma) < tol_ ) return true;

  // move r1p along e_ by -sigma
  r1p -=  sigma * e_;
  pr1p[0] = r1p.x;
  pr1p[1] = r1p.y;
  pr1p[2] = r1p.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool PlaneConstraint::enforce_v(const vector<vector<double> > &r0,
vector<vector<double> > &v0) const
{
  double* pv1 = &v0[is1_][3*ia1_];
  D3vector v1(pv1);

  // vperp: component of v1 perpendicular to the plane
  double vperp = v1 * e_;

  if ( fabs(vperp) < tol_ ) return true;

  v1 -= vperp * e_;

  pv1[0] = v1.x;
  pv1[1] = v1.y;
  pv1[2] = v1.z;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
void PlaneConstraint::compute_force(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  const double* pr1 = &r0[is1_][3*ia1_];
  D3vector r1(pr1);
  const double* pf1 = &f[is1_][3*ia1_];
  D3vector f1(pf1);

  // A positive force on the constraint tends to increase the distance
  force_ = f1 * e_;
}

////////////////////////////////////////////////////////////////////////////////
ostream& PlaneConstraint::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << " <constraint name=\"" << name();
  os << "\" type=\"" << type() << "\">";
  os << " <atoms> " << name1_ << " </atoms>\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "  <velocity> " << setprecision(8) << velocity_ << " </velocity>";
  os << " <weight> " << setprecision(8) << weight_ << " </weight>\n";
  os << "  <value> " << setprecision(8) << distance_ << " </value>";
  os << " <force> " << setprecision(8) << force_ << " </force>\n";
  os << " </constraint>";
  return os;
}
