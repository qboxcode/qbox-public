////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.2 2004-03-11 21:52:31 fgygi Exp $

#ifndef SDIONICSTEPPER_H
#define SDIONICSTEPPER_H

#include "IonicStepper.h"

class SDIonicStepper : public IonicStepper
{
  private:
  
  public:
  
  SDIonicStepper(Sample& s) : IonicStepper(s) {}

  void compute_rp(const vector<vector< double> >& f0);
  void compute_v0(const vector<vector< double> >& f0) {}
  void update_r(void);
  void update_v(void);
  double ekin(void) const { return 0.0; }
  double temp(void) const { return 0.0; } 
};

#endif
