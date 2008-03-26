////////////////////////////////////////////////////////////////////////////////
//
// CGIonicStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CGIonicStepper.C,v 1.1 2008-03-26 04:57:54 fgygi Exp $

#include "CGIonicStepper.h"
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void CGIonicStepper::compute_r(double e0, const vector<vector< double> >& f0)
{
  bool wolfe1;
  if ( first_step_ )
  {
     wolfe1 = true;
     alpha_ = 1.0;
  }
  else
  {
    // Check if the first Wolfe condition is satisfied
    // compute predicted decrease for previous step
    // pred = pc_ * ( r0 - rc_ )
    double pred = 0.0;
    for ( int is = 0; is < r0_.size(); is++ )
    {
      for ( int i = 0; i < r0_[is].size(); i++ )
      {
        pred += pc_[is][i] * ( r0_[is][i] - rc_[is][i] );
      }
    }
    const double sigma_wolfe = 0.1;
    wolfe1 = ( e0 < ec_ - sigma_wolfe * pred );
    //cout << "CGIonicStepper: pred = " << pred << endl;
    //cout << "CGIonicStepper: ec: " << ec_ << endl;
    //cout << "CGIonicStepper: required energy: " << ec_-sigma_wolfe * pred
    //     << " actual: " << e0 << endl;
    //cout << "CGIonicStepper: wolfe1 = " << wolfe1 << endl;
  }

  if ( wolfe1 )
  {
    // actual decrease from last step is sufficient
    // accept r0 as new current point
    rc_ = r0_;
    ec_ = e0;

    // define new descent direction
    if ( first_step_ )
    {
      pc_ = f0;
    }
    else
    {
      // Polak-Ribiere definition
      double num = 0.0, den = 0.0;
      for ( int is = 0; is < r0_.size(); is++ )
      {
        for ( int i = 0; i < r0_[is].size(); i++ )
        {
          const double fctmp = fc_[is][i];
          const double f0tmp = f0[is][i];
          num += f0tmp * ( f0tmp - fctmp );
          den += fctmp * fctmp;
        }
      }
      const double beta = den > 0.0 ? num/den : 0.0;
      //cout << "CGIonicStepper: beta = " << beta << endl;
      for ( int is = 0; is < r0_.size(); is++ )
      {
        for ( int i = 0; i < r0_[is].size(); i++ )
        {
          pc_[is][i] = beta * pc_[is][i] + f0[is][i];
        }
      }
    }
    fc_ = f0;
    alpha_ *= 1.05;
  }
  else
  {
     // backtrack
     alpha_ *= 0.8;
  }
  //cout << "CGIonicStepper: alpha = " << alpha_ << endl;
  for ( int is = 0; is < r0_.size(); is++ )
  {
    for ( int i = 0; i < r0_[is].size(); i++ )
    {
      rp_[is][i] = rc_[is][i] + alpha_ * pc_[is][i];
    }
  }
  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);

  first_step_ = false;
}
