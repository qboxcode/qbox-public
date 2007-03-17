////////////////////////////////////////////////////////////////////////////////
//
// SDCellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDCellStepper.h,v 1.2 2007-03-17 01:14:00 fgygi Exp $

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
