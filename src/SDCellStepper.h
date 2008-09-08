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
// SDCellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDCellStepper.h,v 1.5 2008-09-08 15:56:19 fgygi Exp $

#ifndef SDCELLSTEPPER_H
#define SDCELLSTEPPER_H

#include "CellStepper.h"

class SDCellStepper : public CellStepper
{
  private:

  public:

  SDCellStepper(Sample& s) : CellStepper(s) {}

  void compute_new_cell(const std::valarray<double>& sigma_eks);
  void update_cell(void);
  double ekin(void) const { return 0.0; }
};
#endif
