////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDIonicStepper.C,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#include "SDIonicStepper.h"

void SDIonicStepper::update(const vector<vector< double> >& fion)
{
  // Steepest descent step
  for ( int is = 0; is < tau0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < tau0_[is].size(); i++ )
    {
        tau0_[is][i] += dt2bym * fion[is][i];
    }
  }
  atoms_.set_positions(tau0_);
}
