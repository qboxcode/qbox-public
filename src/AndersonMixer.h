////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.h,v 1.2 2004-12-10 01:04:06 fgygi Exp $

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <valarray>
#include <cassert>
using namespace std;

#include "Context.h"

class AndersonMixer
{
  int     n_;                    // size of vectors
  const   Context* const pctxt_; // pointer to relevant Context, null if local
  double  theta_max_;
 
  valarray<double> flast_;       // last residual
  bool extrapolate_;             // state variable

  public:

  AndersonMixer(const int n, const Context* const pctxt) :
    n_(n), pctxt_(pctxt), extrapolate_(false), theta_max_(2.0)
  {
    assert( n > 0 );
    flast_.resize(n);
  }

  void update(const double* f, double* theta, double* fbar);
  void restart(void);
  void set_theta_max(double theta_max) { theta_max_ = theta_max; }
  double theta_max(void) const { return theta_max_; }
};
#endif
