////////////////////////////////////////////////////////////////////////////////
//
// CellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellStepper.h,v 1.2 2007-03-17 01:14:00 fgygi Exp $

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
