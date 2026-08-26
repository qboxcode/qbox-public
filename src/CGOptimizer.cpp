////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// CGOptimizer.cpp
//
////////////////////////////////////////////////////////////////////////////////
#include "CGOptimizer.h"
#include <iostream>
#include <cassert>
#include <algorithm>
using namespace std;
////////////////////////////////////////////////////////////////////////////////
void CGOptimizer::compute_xp(const valarray<double>& x, const double f,
                             valarray<double>& g, valarray<double>& xp)
{
  // Use the function value f and the gradient g at x to generate a new point xp
  // using the Polak-Ribiere+ CG algorithm
  // return xp=x if the 2-norm of g is smaller than tol
  const double tol = 1.0e-18;

  assert(x.size()==n_ && g.size()==n_ && xp.size()==n_);

  double fp;
  // define the descent direction
  if ( first_step_ )
  {
    p_ = -g;
    gm_ = g;

    x0_ = x;
    f0_ = f;

    g0norm2_ = (g*g).sum();
    if ( g0norm2_ < tol )
    {
      xp = x;
      return;
    }
    fp = -g0norm2_;
    fp0_ = fp;
    linmin_.reset();
    // The following call initializes linmin and returns alpha_start
    alpha_ = linmin_.next_alpha(0.0,f,fp);
    if ( debug_print )
      cout << "  CGOptimizer: first_step: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;

    xp = x0_ + alpha_ * p_;
    first_step_ = false;
    return;
  }

  // This is not the first CG step
  // fp: derivative along the current descent direction p_
  // fp = df(x0+alpha*p)/dalpha at x
  fp = (g*p_).sum();
  alpha_ = linmin_.next_alpha(alpha_,f,fp);
  if ( debug_print )
    cout << "  CGOptimizer: alpha=" << alpha_
         << " f=" << f << " fp=" << fp << endl;

  if ( linmin_.fail() )
  {
    // line minimization failed
    if ( debug_print )
      cout << "  CGOptimizer: line minimization failed" << endl;

    // restart from current point
    p_ = -g;
    gm_ = g;

    x0_ = x;
    f0_ = f;

    g0norm2_ = (g*g).sum();
    if ( g0norm2_ < tol )
    {
      xp = x;
      return;
    }
    fp = -g0norm2_;
    fp0_ = fp;

    linmin_.reset();
    // The following call initializes linmin and returns alpha_start
    alpha_ = linmin_.next_alpha(0.0,f,fp);

    if ( debug_print )
      cout << "  CGOptimizer: restart after fail: alpha=" << alpha_
         << " f=" << f << " fp=" << fp << endl;

    xp = x0_ + alpha_ * p_;
    first_step_ = false;
    return;
  }

  if ( linmin_.reached_alpha_max() )
  {
    if ( debug_print )
      cout << "  CGOptimizer: linmin reached alpha_max" << endl;

    // continue search in the same direction

    x0_ = x;
    f0_ = f;
    fp = (g*p_).sum();
    g0norm2_ = (g*g).sum();
    gm_ = g;
    fp0_ = fp;

    // choose value of alpha for the next step
    alpha_ = alpha_max_ratio_ * linmin_.alpha_max();

    // restart linmin with alpha_
    linmin_.reset();
    // the following call initializes linmin and returns alpha_
    alpha_ = linmin_.next_alpha(alpha_,f,fp);

    if ( debug_print )
      cout << "  CGOptimizer: restart: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;
    xp = x0_ + alpha_ * p_;
    return;
  }

  if ( linmin_.done() )
  {
    // wolfe1_ && wolfe2_ are true at alpha_
    if ( debug_print )
      cout << "  CGOptimizer: done with current descent direction" << endl;
    // define a new descent direction p_ using the Polak-Ribiere+ formula
    assert(g0norm2_ > 0.0);

    // Polak-Ribiere+: clamp beta >= 0 to ensure descent direction
    double beta = ((g*g).sum()-(gm_*g).sum()) / g0norm2_;
    beta = max(beta, 0.0);

    // clamp above at beta_max if set
    if ( beta_max_ > 0.0 && beta > beta_max_ )
      beta = beta_max_;

    if ( debug_print )
      cout << "  CGOptimizer: beta = " << beta << endl;
    p_ = beta * p_ - g;

    x0_ = x;
    f0_ = f;
    fp = (g*p_).sum();
    g0norm2_ = (g*g).sum();
    gm_ = g;
    fp0_ = fp;

    // estimate of value of alpha for the next line search
    alpha_ = min(alpha_, alpha_max_ratio_ * linmin_.alpha_max());

    // restart linmin with alpha_
    linmin_.reset();
    // the following call initializes linmin and returns alpha_
    alpha_ = linmin_.next_alpha(alpha_,f,fp);

    if ( debug_print )
      cout << "  CGOptimizer: restart: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;
    xp = x0_ + alpha_ * p_;
    return;
  }

  // normal case
  xp = x0_ + alpha_ * p_;
}
