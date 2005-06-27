////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.4 2005-06-27 22:20:10 fgygi Exp $

#ifndef SDIONICSTEPPER_H
#define SDIONICSTEPPER_H

#include "IonicStepper.h"

class SDIonicStepper : public IonicStepper
{
  private:
  
  public:
  
  SDIonicStepper(Sample& s) : IonicStepper(s) {}

  void compute_r(double e0, const vector<vector< double> >& f0);
  void compute_v(double e0, const vector<vector< double> >& f0) {}
};

#endif
