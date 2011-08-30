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

#ifndef SDCELLSTEPPER_H
#define SDCELLSTEPPER_H

#include "CellStepper.h"

class SDCellStepper : public CellStepper
{
  private:

  public:

  SDCellStepper(Sample& s);

  void compute_new_cell(double e0, const std::valarray<double>& sigma,
    const std::vector<std::vector< double> >& f0);
  void update_cell(void);
  double ekin(void) const { return 0.0; }
};
#endif
