////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.3 2004-12-17 23:37:25 fgygi Exp $

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
};

#endif
