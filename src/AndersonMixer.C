////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.C,v 1.1 2004-12-02 22:24:16 fgygi Exp $

#include "AndersonMixer.h"
#include "blas.h"

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::restart(void)
{
  extrapolate_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::update(const double* f, double* theta, double* fbar)
{
    *theta = 0.0;
    if ( extrapolate_ )
    {
      int ione=1;
      valarray<double> tmp0(n_);

      // compute theta = - a / b
      // tmp0 = delta_F = f - flast
      for ( int i = 0; i < n_; i++ )
        tmp0[i] = f[i] - flast_[i];
      
      // a = delta_F * F
      double a = ddot(&n_, &tmp0[0], &ione, f, &ione);
        
      // b = delta_F * delta_F
      double b = ddot(&n_, &tmp0[0], &ione, &tmp0[0], &ione);
      
      double work[2] = { a, b };
      ctxt_.dsum(2,1,work,2);
      a = work[0];
      b = work[1];

      if ( b != 0.0 )
        *theta = - a / b;
      else
        *theta = 0.0;
      
      // test if residual has increased
      if ( *theta <= -1.0 )
      {
        *theta = 1.0;
      }
      
      *theta = min(theta_max_,*theta);
    }
    
    // fbar = f + theta * ( f - flast )
    // flast_ = f_
    for ( int i = 0; i < n_; i++ )
    {
      const double ftmp = f[i];
      fbar[i] = ftmp + *theta * ( ftmp - flast_[i] );
      flast_[i] = ftmp;
    } 
    extrapolate_ = true;
}
