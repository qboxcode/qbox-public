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
// $Id: CellStepper.h,v 1.5 2008-09-08 15:56:18 fgygi Exp $

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

  public:

  CellStepper (Sample& s) : s_(s), atoms_(s.atoms), ekin_(0.0) {}

  virtual void compute_new_cell(const std::valarray<double>& sigma) = 0;
  virtual void update_cell(void) = 0;

  double ekin(void) const { return ekin_; }
  virtual ~CellStepper() {}
};
#endif
