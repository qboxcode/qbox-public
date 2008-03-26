////////////////////////////////////////////////////////////////////////////////
//
// CGIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CGIonicStepper.h,v 1.1 2008-03-26 04:57:54 fgygi Exp $

#ifndef CGIONICSTEPPER_H
#define CGIONICSTEPPER_H

#include "IonicStepper.h"
#include <vector>

class CGIonicStepper : public IonicStepper
{
  private:

  bool first_step_;
  std::vector<std::vector< double> > rc_;
  std::vector<std::vector< double> > pc_;
  std::vector<std::vector< double> > fc_;
  double ec_;
  double alpha_;

  public:

  CGIonicStepper(Sample& s) : IonicStepper(s), first_step_(true) {}

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
};

#endif
