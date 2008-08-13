////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License 
// as published by the Free Software Foundation, either version 2 of 
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.h,v 1.6 2008-08-13 06:39:42 fgygi Exp $

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
