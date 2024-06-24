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
// PlaneConstraint.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PLANECONSTRAINT_H
#define PLANECONSTRAINT_H

#include "Constraint.h"
#include "D3vector.h"
#include <cassert>
#include <cmath> // fabs

class AtomSet;

class PlaneConstraint : public Constraint
{
  std::string name1_;
  int ia1_, is1_;
  double distance_, velocity_, force_, weight_, tol_;
  D3vector e_;

  public:

  PlaneConstraint(std::string name, std::string name1, D3vector p,
                  double distance, double velocity, double tolerance):
  name1_(name1), distance_(distance), e_(normalized(p)),
  velocity_(velocity), tol_(tolerance)
  {
    name_ = name;
    names_.resize(1);
    names_[0] = name1_;
    force_ = 0.0;
    weight_ = 1.0;
  }
  ~PlaneConstraint(void) {}

  std::string type(void) const { return "plane"; }
  double value(void) const { return distance_; }
  double velocity(void) const { return velocity_; }
  double force(void) const { return force_; }
  double weight(void) const { return weight_; }
  double tolerance(void) const { return tol_; }
  void set_value(double value)
  {
    distance_ = value;
  }
  void set_velocity(double velocity)
  {
    velocity_ = velocity;
  }

  void setup(const AtomSet& atoms);
  int dofs(void) const { return 1; }
  void update(double dt);
  bool enforce_r(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &rp) const;
  bool enforce_v(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &v0) const;
  void compute_force(const std::vector<std::vector<double> > &r0,
                     const std::vector<std::vector<double> > &f);
  std::ostream& print( std::ostream& os );

};
#endif
