////////////////////////////////////////////////////////////////////////////////
//
// MDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIonicStepper.C,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#include "MDIonicStepper.h"

////////////////////////////////////////////////////////////////////////////////
double MDIonicStepper::temp(void) const
{
  const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
  if ( ndofs_ > 0.0 )
    return 2.0 * ( ekin_ / boltz ) / ndofs_;
  else
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::update(const vector<vector< double> >& fion)
{
  // Verlet update
  if ( thermostat_ )
  {
    eta_ = 0.0;
    double ekin_ion = 0.0;
    // compute damping factor eta
    // compute ekin_ion and temp_ion before step using a first order
    // approximation. Note: ekin_ is recomputed after the step using
    // a second-order approximation.
    for ( int is = 0; is < tau0_.size(); is++ )
    {
      const double dt2bym = dt_ * dt_ / pmass_[is];
      if ( dt_ != 0.0 )
      {
        for ( int i = 0; i < tau0_[is].size(); i++ )
        {
          const double v = ( tau0_[is][i] - taum_[is][i] ) / dt_;
          ekin_ion += 0.5 * pmass_[is] * v * v;
        }
      }
    }
    // Next line: linear thermostat, width of tanh = 100.0
    eta_ = tanh ( ( temp() - th_temp_ ) / 100. ) / th_time_;
  }
 
  // make Verlet step with damping eta, and recompute ekin_ to second order
  ekin_ = 0.0;
  for ( int is = 0; is < tau0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < tau0_[is].size(); i++ )
    {
      const double taup = tau0_[is][i] +
      (1.0 - dt_*eta_) * (tau0_[is][i] - taum_[is][i]) +
      dt2bym * fion[is][i];
      if ( dt_ != 0.0 )
      {
        const double v = 0.5 * ( taup - taum_[is][i] ) / dt_;
        ekin_ += 0.5 * pmass_[is] * v * v;
        vel_[is][i] = v;
      }
      taum_[is][i] = tau0_[is][i];
      tau0_[is][i] = taup;
    }
  }
  
  atoms_.set_positions(tau0_);
  atoms_.set_velocities(vel_);
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::stoermer_start(const vector<vector< double> >& fion)
{
  // compute taum from tau0,vel,fion and dt before first step
  // First step of Stoermer's rule
  // x1 = x0 + h * v0 + 0.5 * dt2/m * f0
  
  // x-1 = x0 - h * v0 + 0.5 * dt2/m * f0
 
  atoms_.get_velocities(vel_);
 
  for ( int is = 0; is < tau0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < tau0_[is].size(); i++ )
    {
      const double taum = tau0_[is][i] - dt_ * vel_[is][i] +
                          0.5 * dt2bym * fion[is][i];

      taum_[is][i] = taum;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void MDIonicStepper::stoermer_end(const vector<vector< double> >& fion)
{
  // compute vel from tau0, taum, fion and dt at last step
  // Last step of Stoermer's rule
  // v = (xn - x(n-1) )/dt + 0.5 * dt/m * fn
 
  if ( dt_ != 0.0 )
  {
    for ( int is = 0; is < tau0_.size(); is++ )
    {
      const double dtbym = dt_ / pmass_[is];
      for ( int i = 0; i < tau0_[is].size(); i++ )
      {
        vel_[is][i] = ( tau0_[is][i] - taum_[is][i] ) / dt_ +
                     0.5 * dtbym * fion[is][i];
      }
    }
  }
  atoms_.set_velocities(vel_);
}
