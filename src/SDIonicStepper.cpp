////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SDIonicStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SDIonicStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  // Steepest descent step
  for ( int is = 0; is < r0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] + dt2bym * f0[is][i];
    }
  }

  if ( s_.ctrl.lock_cm )
    reset_rcm(r0_,rp_);

  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
}
