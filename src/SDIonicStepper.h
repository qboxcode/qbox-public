////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#ifndef SDIONICSTEPPER_H
#define SDIONICSTEPPER_H

#include "IonicStepper.h"

class SDIonicStepper : public IonicStepper
{
  private:
  
  public:
  
  SDIonicStepper(Sample& s) : IonicStepper(s) {}

  void update(const vector<vector< double> >& fion);
  double ekin(void) const { return 0.0; }
  double temp(void) const { return 0.0; } 
};

#endif
