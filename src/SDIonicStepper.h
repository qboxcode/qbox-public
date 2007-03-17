////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.5 2007-03-17 01:14:00 fgygi Exp $

#ifndef SDIONICSTEPPER_H
#define SDIONICSTEPPER_H

#include "IonicStepper.h"
#include <vector>

class SDIonicStepper : public IonicStepper
{
  private:
  
  public:
  
  SDIonicStepper(Sample& s) : IonicStepper(s) {}

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};

#endif
