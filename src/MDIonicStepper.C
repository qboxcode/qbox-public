////////////////////////////////////////////////////////////////////////////////
//
// MDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIonicStepper.C,v 1.7 2004-05-20 00:20:26 fgygi Exp $

#include "MDIonicStepper.h"

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::compute_rp(const vector<vector< double> >& f0)
{
  // f0 contains forces at r0
  // Compute new positions rp using the Verlet algorithm
  // compute velocity at t+1/2 vhalf_
  
  ekin_ = 0.0;
  atoms_.get_velocities(v0_);
  for ( int is = 0; is < v0_.size(); is++ )
  {
    const double dtby2m = dt_ / ( 2.0 * pmass_[is] );
    if ( dt_ != 0.0 )
    {
      for ( int i = 0; i < v0_[is].size(); i++ )
      {
        const double v = v0_[is][i] + dtby2m * f0[is][i];
        ekin_ += 0.5 * pmass_[is] * v0_[is][i] * v0_[is][i];
        vhalf_[is][i] = v;
      }
    }
  }
  
  // ekin_ is the kinetic energy computed from v0
  
  if ( thermostat_ )
  {
    eta_ = tanh ( ( temp() - th_temp_ ) / th_width_ ) / th_time_;
    if ( s_.ctxt_.onpe0() )
    {
      cout << "  <!-- Thermostat: temp=" << temp() << " -->" << endl;
      cout << "  <!-- Thermostat: tref=" << th_temp_ << " -->" << endl;
      cout << "  <!-- Thermostat: eta=" << eta_ << " -->" << endl;
    }
    
    const double fac = (1.0 - eta_ * dt_);
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        vhalf_[is][i] *= fac;
      }
    }
  }
 
  // compute rp
  atoms_.get_positions(r0_);
  for ( int is = 0; is < r0_.size(); is++ )
  {
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] + vhalf_[is][i] * dt_;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::compute_v0(const vector<vector< double> >& f0)
{
  // compute velocities v0_ using vhalf_ and force f0(r0)
  // Note: vhalf contains the velocity at t-1/2 since compute_v0 is called
  // after a call to update_r()
  
  ekin_ = 0.0;
  for ( int is = 0; is < vhalf_.size(); is++ )
  {
    assert(pmass_[is] > 0.0);
    const double dtby2m = dt_ / ( 2.0 * pmass_[is] );
    for ( int i = 0; i < vhalf_[is].size(); i++ )
    {
      const double v = vhalf_[is][i] + dtby2m * f0[is][i];
      ekin_ += 0.5 * pmass_[is] * v * v;
      v0_[is][i] = v;
    }
  }
  // ekin_ contains the kinetic energy computed from v0
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::update_r(void)
{
  r0_ = rp_;
  atoms_.set_positions(r0_);
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::update_v(void)
{
  atoms_.set_velocities(v0_);
}
