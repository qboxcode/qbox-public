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
// BMDIonicStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////
#include "BMDIonicStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  // f0 contains forces at r0
  // e0 contains energy at r0
  // Compute new positions rp using the velocity Verlet algorithm
  // enforce constraints for rp
  // update rm <- r0, r0 <- rp, and update atomset

  // compute rp
  for ( int is = 0; is < r0_.size(); is++ )
    for ( int i = 0; i < r0_[is].size(); i++ )
      rp_[is][i] = r0_[is][i] + v0_[is][i] + bmd_fac_ * f0[is][i];

  if ( s_.ctrl.lock_cm )
    reset_rcm(r0_,rp_);

  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  fm_ = f0;
  em_ = e0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
}

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_v(double e0, const vector<vector< double> >& f0)
{
  // compute velocities v0_ using r0_, rm_ and f0(r0)
  // enforce constraints for vp
  // adjust velocities with the thermostat

  e0_ = e0;
  for ( int is = 0; is < v0_.size(); is++ )
    for ( int i = 0; i < v0_[is].size(); i++ )
      v0_[is][i] = r0_[is][i] - rm_[is][i] + bmd_fac_ * f0[is][i];

  // check if energy increased
  if ( e0_ > em_ )
  {
    // rescale velocity or reset to zero
    for ( int is = 0; is < r0_.size(); is++ )
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        // stop if velocity component is opposed to force
        if ( v0_[is][i] * f0[is][i] < 0.0 )
          v0_[is][i] = 0.0;
      }
  }
  else
  {
    // accelerate
    for ( int is = 0; is < r0_.size(); is++ )
      for ( int i = 0; i < r0_[is].size(); i++ )
        v0_[is][i] *= 1.05;
  }

  if ( s_.ctrl.lock_cm )
    reset_vcm(v0_);

  constraints_.enforce_v(r0_,v0_);
  compute_ekin();
}

////////////////////////////////////////////////////////////////////////////////
void BMDIonicStepper::compute_ekin(void)
{
  ekin_ = 0.0;
  for ( int is = 0; is < v0_.size(); is++ )
  {
    for ( int i = 0; i < v0_[is].size(); i++ )
    {
      const double v = v0_[is][i];
      ekin_ += v * v / ( 4.0 * bmd_fac_ );
    }
  }
}
