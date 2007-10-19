////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.h,v 1.5 2007-10-19 16:24:03 fgygi Exp $

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <valarray>
#include <cassert>
#include "Context.h"

class AndersonMixer
{
  int     n_;                    // size of vectors
  const   Context* const pctxt_; // pointer to relevant Context, null if local
  double  theta_max_; // maximum extrapolation
  double  theta_nc_;  // negative curvature value

  std::valarray<double> flast_;       // last residual
  bool extrapolate_;             // state variable

  public:

  AndersonMixer(const int n, const Context* const pctxt) :
    n_(n), pctxt_(pctxt), extrapolate_(false), theta_max_(2.0), theta_nc_(0.0)
  {
    assert( n > 0 );
    flast_.resize(n);
  }

  void update(const double* f, double* theta, double* fbar);
  void restart(void);
  void set_theta_max(double theta_max) { theta_max_ = theta_max; }
  void set_theta_nc(double theta_nc) { theta_nc_ = theta_nc; }
  double theta_max(void) const { return theta_max_; }
  double theta_nc(void) const { return theta_nc_; }
};
#endif
