////////////////////////////////////////////////////////////////////////////////
//
// SDAIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDAIonicStepper.h,v 1.6 2007-03-17 01:14:00 fgygi Exp $

#ifndef SDAIONICSTEPPER_H
#define SDAIONICSTEPPER_H

#include "IonicStepper.h"
#include "AndersonMixer.h"

class SDAIonicStepper : public IonicStepper
{
  private:
  
  std::vector<double> f_;
  std::vector<double> fbar_;
  double theta_;
  bool first_step_;
  AndersonMixer mixer_;
  
  public:
  
  SDAIonicStepper(Sample& s) : IonicStepper(s), first_step_(true), theta_(0),
  mixer_(3*atoms_.size(), 0)
  {
    f_.resize(3*atoms_.size());
    fbar_.resize(3*atoms_.size());
    rm_ = r0_;
    mixer_.set_theta_max(1.0);
    //mixer_.set_theta_nc(-0.5); // default: 0.0
  }

  void compute_r(double e0, const std::vector<std::vector< double> >& f0);
  void compute_v(double e0, const std::vector<std::vector< double> >& f0) {}
  void reset(void) { first_step_ = true; }
};

#endif
