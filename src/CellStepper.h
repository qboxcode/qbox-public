////////////////////////////////////////////////////////////////////////////////
//
// CellStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellStepper.h,v 1.1 2004-03-11 21:58:10 fgygi Exp $

#ifndef CELLSTEPPER_H
#define CELLSTEPPER_H

#include "Sample.h"
#include <valarray>
#include <iostream>
#include <iomanip>
using namespace std;

class CellStepper
{
  protected:
  
  Sample& s_;
  AtomSet& atoms_;
  double ekin_;
  UnitCell cellp;

  public:
  
  CellStepper (Sample& s) : s_(s), atoms_(s.atoms), ekin_(0.0) {}
  
  virtual void compute_new_cell(const valarray<double>& sigma) = 0;
  virtual void update_cell(void) = 0;
  
  double ekin(void) const { return ekin_; }
  virtual ~CellStepper() {}
};
#endif
