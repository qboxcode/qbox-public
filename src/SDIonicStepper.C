////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.C,v 1.2 2004-03-11 21:52:32 fgygi Exp $

#include "SDIonicStepper.h"

////////////////////////////////////////////////////////////////////////////////
void SDIonicStepper::compute_rp(const vector<vector< double> >& f0)
{
  // Steepest descent step
  atoms_.get_positions(r0_);
  for ( int is = 0; is < r0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] + dt2bym * f0[is][i];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void SDIonicStepper::update_r(void)
{
  r0_ = rp_;
  atoms_.set_positions(r0_);
}

////////////////////////////////////////////////////////////////////////////////
void SDIonicStepper::update_v(void)
{
  atoms_.reset_velocities(); // set velocities to zero
}
