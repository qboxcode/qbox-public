////////////////////////////////////////////////////////////////////////////////
//
// SDCellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDCellStepper.h,v 1.1 2004-03-11 21:58:10 fgygi Exp $

#ifndef SDCELLSTEPPER_H
#define SDCELLSTEPPER_H

#include "CellStepper.h"

class SDCellStepper : public CellStepper
{
  private:
  
  public:
  
  SDCellStepper(Sample& s) : CellStepper(s) {}

  void compute_new_cell(const valarray<double>& sigma_eks);
  void update_cell(void);
  double ekin(void) const { return 0.0; }
};
#endif
