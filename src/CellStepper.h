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
// CellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CELLSTEPPER_H
#define CELLSTEPPER_H

#include "Sample.h"
#include <valarray>

class CellStepper
{
  protected:

  Sample& s_;
  AtomSet& atoms_;
  double ekin_;
  UnitCell cellp;
  std::vector<double> u_, up_;

  public:

  CellStepper (Sample& s) : s_(s), atoms_(s.atoms), ekin_(0.0) {}

  virtual void compute_new_cell(double e0,const std::valarray<double>& sigma,
    const std::vector<std::vector< double> >& f0) = 0;
  void enforce_constraints(double* u);
  virtual void update_cell(void) = 0;

  double ekin(void) const { return ekin_; }
  virtual ~CellStepper() {}
};
#endif
