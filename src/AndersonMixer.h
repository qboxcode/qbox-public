////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.h,v 1.1 2004-12-02 22:24:16 fgygi Exp $

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <valarray>
#include <cassert>
using namespace std;

#include "Context.h"

class AndersonMixer
{
  int     n_;              // size of vectors
  const   Context ctxt_;
  double  theta_max_;
 
  valarray<double> flast_; // last residual
  bool extrapolate_;       // state variable

  public:
    
  AndersonMixer(const int n, const Context& ctxt) :
    n_(n), ctxt_(ctxt), extrapolate_(false), theta_max_(2.0)
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
