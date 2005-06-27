////////////////////////////////////////////////////////////////////////////////
//
// MDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIonicStepper.C,v 1.9 2005-06-27 22:28:33 fgygi Exp $

#include "MDIonicStepper.h"

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  // f0 contains forces at r0
  // Compute new positions rp using the velocity Verlet algorithm
  // enforce constraints for rp
  // update rm <- r0, r0 <- rp, and update atomset
  
  // compute rp
  for ( int is = 0; is < r0_.size(); is++ )
  {
    const double dt2by2m = dt_ * dt_ / ( 2.0 * pmass_[is] );
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] + v0_[is][i] * dt_ + dt2by2m * f0[is][i];
    }
  }
  
  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::compute_v(double e0, const vector<vector< double> >& f0)
{
  // compute velocities v0_ using r0_, rm_ and f0(r0)
  // enforce constraints for vp
  // adjust velocities with the thermostat
  
  assert(dt_ > 0.0);
  for ( int is = 0; is < v0_.size(); is++ )
  {
    assert(pmass_[is] > 0.0);
    const double dtby2m = dt_ / ( 2.0 * pmass_[is] );
    for ( int i = 0; i < v0_[is].size(); i++ )
    {
      const double vhalf = ( r0_[is][i] - rm_[is][i] ) / dt_;
      v0_[is][i] = vhalf + dtby2m * f0[is][i];
    }
  }
  compute_ekin();
  if ( thermostat_ )
  {
    eta_ = tanh ( ( temp() - th_temp_ ) / th_width_ ) / th_time_;
    if ( s_.ctxt_.onpe0() )
    {
      cout << "  <!-- thermostat: temp=" << temp() << " -->" << endl;
      cout << "  <!-- thermostat: tref=" << th_temp_ << " -->" << endl;
      cout << "  <!-- thermostat: eta=" << eta_ << " -->" << endl;
    }
    
    const double fac = (1.0 - eta_ * dt_);
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        v0_[is][i] *= fac;
      }
    }
    ekin_ *= fac * fac;
  }
  constraints_.enforce_v(r0_,v0_);
  atoms_.set_velocities(v0_);
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::compute_ekin(void)
{
  ekin_ = 0.0;
  for ( int is = 0; is < v0_.size(); is++ )
  {
    for ( int i = 0; i < v0_[is].size(); i++ )
    {
      const double v = v0_[is][i];
      ekin_ += 0.5 * pmass_[is] * v * v;
    }
  }
}
