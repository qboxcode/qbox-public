////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.h,v 1.7 2008-03-26 04:57:54 fgygi Exp $

#ifndef SDIONICSTEPPER_H
#define SDIONICSTEPPER_H

#include "IonicStepper.h"
#include <vector>

class SDIonicStepper : public IonicStepper
{
  private:

  bool first_step_;
  std::vector<std::vector< double> > rc_;
  std::vector<std::vector< double> > pc_;
  double ec_;
  double alpha_;

  public:

  SDIonicStepper(Sample& s) : IonicStepper(s), first_step_(true) {}

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};

#endif
