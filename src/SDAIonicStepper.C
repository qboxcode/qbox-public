////////////////////////////////////////////////////////////////////////////////
//
// SDAIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SDAIonicStepper.C,v 1.3 2004-12-18 23:21:42 fgygi Exp $

#include "SDAIonicStepper.h"
#include <iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDAIonicStepper::compute_rp(const vector<vector< double> >& f0)
{
  // Steepest descent step
  atoms_.get_positions(r0_);
  
  if ( first_step_ )
  {
    rm_ = r0_;
  }
  
  // copy forces f0 in linear array f_
  int k = 0;
  double sum = 0.0;
  for ( int is = 0; is < f0.size(); is++ )
  {
    for ( int i = 0; i < f0[is].size(); i++ )
    {
      f_[k++] = f0[is][i];
      sum += f0[is][i]*f0[is][i];
    }
  }
  if ( s_.ctxt_.onpe0() )
    cout << "<sda_residual> " << sum 
         << "</sda_residual>" << endl;
  
  mixer_.update(&f_[0],&theta_,&fbar_[0]);
  if ( s_.ctxt_.onpe0() )
    cout << "<sda_theta> " << theta_ 
         << "</sda_theta>" << endl;
    
  k = 0;
  for ( int is = 0; is < r0_.size(); is++ )
  {
    const double dt2bym = dt_ * dt_ / pmass_[is];
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = r0_[is][i] +
                   theta_ * ( r0_[is][i] - rm_[is][i] ) +
                   dt2bym * fbar_[k++];
    }
  }
  first_step_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void SDAIonicStepper::update_r(void)
{
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
}

////////////////////////////////////////////////////////////////////////////////
void SDAIonicStepper::update_v(void)
{
  atoms_.reset_velocities(); // set velocities to zero
}
