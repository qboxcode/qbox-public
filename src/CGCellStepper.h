////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// CGCellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CGCELLSTEPPER_H
#define CGCELLSTEPPER_H

#include "CellStepper.h"
#include "LineMinimizer.h"

class CGCellStepper : public CellStepper
{
  private:

  bool first_step_;
  int nat_;
  std::valarray<double> xc_, pc_, fc_;
  std::vector<std::vector<double> > rc_, rp_;
  double ec_, fpc_;
  double alpha_, sigma1_, sigma2_;
  LineMinimizer linmin_;

  public:

  CGCellStepper(Sample& s);
  void compute_new_cell(double e0, const std::valarray<double>& sigma,
                        const std::vector<std::vector< double> >& fion);
  void update_cell(void);
  double ekin(void) const { return 0.0; }
};
#endif
