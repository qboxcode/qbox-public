////////////////////////////////////////////////////////////////////////////////
//
// SDAIonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDAIonicStepper.h,v 1.1 2004-12-10 01:13:04 fgygi Exp $

#ifndef SDAIONICSTEPPER_H
#define SDAIONICSTEPPER_H

#include "IonicStepper.h"
#include "AndersonMixer.h"

class SDAIonicStepper : public IonicStepper
{
  private:
  
  vector<vector<double> > rm_;
  vector<double> f_;
  vector<double> fbar_;
  double theta_;
  bool first_step_;
  AndersonMixer mixer_;
  
  public:
  
  SDAIonicStepper(Sample& s) : IonicStepper(s), first_step_(true), theta_(0),
  mixer_(ndofs_, 0)
  {
    f_.resize(ndofs_);
    fbar_.resize(ndofs_);
    rm_ = r0_;
  }

  void compute_rp(const vector<vector< double> >& f0);
  void compute_v0(const vector<vector< double> >& f0) {}
  void update_r(void);
  void update_v(void);
  double ekin(void) const { return 0.0; }
  double temp(void) const { return 0.0; }
  void reset(void) { first_step_ = true; }
};

#endif
