////////////////////////////////////////////////////////////////////////////////
//
// MDIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MDIonicStepper.C,v 1.11 2006-11-05 02:06:29 fgygi Exp $

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
  
  if ( thermostat_ == "SCALING" )
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
  else if ( thermostat_ == "ANDERSEN" )
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    for ( int is = 0; is < nsp_; is++ )
    {
      const double m = pmass_[is];        
      const double width = sqrt( boltz * th_temp_ / m );
      for ( int ia = 0; ia < na_[is]; ia++ )  
      {                                                 
        // if th_time is zero, set probability to one
        if ( th_time_ == 0.0 || drand48() < dt_ / th_time_ )               
        {                                               
          // draw gaussian random variables                        
          double xi0 = drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() - 6.0;    
          double xi1 = drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() - 6.0;    
          double xi2 = drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() +         
                       drand48() + drand48() + drand48() - 6.0;    
          v0_[is][3*ia+0] = width * xi0;    
          v0_[is][3*ia+1] = width * xi1;    
          v0_[is][3*ia+2] = width * xi2;    
        }
      }
    }
  }
  else if ( thermostat_ == "LOWE" )
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    
    // scan all atom pairs in the space (is,ia)
    //int npairs = 0;
    for ( int is1 = 0; is1 < nsp_; is1++ )
    {
      for ( int is2 = is1; is2 < nsp_; is2++ )
      {
        const double m1 = pmass_[is1];        
        const double m2 = pmass_[is2];        
        const double mu = m1*m2/(m1+m2);
        const double width = sqrt(boltz * th_temp_ / mu);
        for ( int ia1 = 0; ia1 < na_[is1]; ia1++ )
        {
          int ia2min = 0;
          if ( is1 == is2 ) ia2min = ia1 + 1;
          for ( int ia2 = ia2min; ia2 < na_[is2]; ia2++ )
          {
            // if th_time is zero, set probability to one
            if ( th_time_ == 0.0 || drand48() < dt_ / th_time_ )               
            {                                               
              //cout << " pair " << is1 << " " << ia1 << " "  
              //     << is2 << " " << ia2 << endl;
              D3vector r1(&r0_[is1][3*ia1]);
              D3vector r2(&r0_[is2][3*ia2]);
              D3vector r12 = r1 - r2;
              s_.wf.cell().fold_in_ws(r12);
              D3vector e12 = normalized(r12);
              D3vector v1(&v0_[is1][3*ia1]);
              D3vector v2(&v0_[is2][3*ia2]);
              D3vector v12 = v1 - v2;
              // draw a gaussian random variable
              double xi = drand48() + drand48() + drand48() +
                          drand48() + drand48() + drand48() +
                          drand48() + drand48() + drand48() +
                          drand48() + drand48() + drand48() - 6.0;
              double lambda = xi * width;
              D3vector dv12 = mu * ( lambda - v12*e12 ) * e12;
              D3vector v1p = v1 + (1.0/m1) * dv12;
              D3vector v2p = v2 - (1.0/m2) * dv12;
              
              v0_[is1][3*ia1+0] = v1p.x;
              v0_[is1][3*ia1+1] = v1p.y;
              v0_[is1][3*ia1+2] = v1p.z;
              
              v0_[is2][3*ia2+0] = v2p.x;
              v0_[is2][3*ia2+1] = v2p.y;
              v0_[is2][3*ia2+2] = v2p.z;             
            }
            //npairs++;                                    
          }
        }
      }
    }
    //cout << " npairs: " << npairs << endl;
    compute_ekin();
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
