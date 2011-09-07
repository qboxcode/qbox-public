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
// CGOptimizer.C
//
////////////////////////////////////////////////////////////////////////////////
#include "CGOptimizer.h"
#include <iostream>
#include <cassert>
#include <numeric>
using namespace std;
////////////////////////////////////////////////////////////////////////////////
double norm2(const std::valarray<double>&v)
{
  double sum = 0.0;
  for ( int i = 0; i < v.size(); i++ )
    sum += v[i]*v[i];
  return sum;
}
////////////////////////////////////////////////////////////////////////////////
void CGOptimizer::compute_xp(const valarray<double>& x, const double f,
                             const valarray<double>& g, valarray<double>& xp)
{
  // Use the function value f and the gradient g at x to generate a new point xp
  // using the Fletcher-Powell CG algorithm
  // return xp=x if the 2-norm of g is smaller than tol
  const double tol = 1.0e-18;

  assert(x.size()==n_ && g.size()==n_ && xp.size()==n_);

  double fp;
  // define the descent direction
  if ( first_step_ )
  {
    p_ = -g;

    x0_ = x;
    f0_ = f;

    g0norm2_ = norm2(g);
    if ( g0norm2_ < tol )
    {
      xp = x;
      return;
    }
    fp = -g0norm2_;
    fp0_ = fp;
    linmin_.reset();
    alpha_ = linmin_.next_alpha(alpha_,f,fp);
    if ( debug_print )
      cout << "  CGOptimizer: first_step: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;

    xp = x0_ + alpha_ * p_;
    first_step_ = false;
  }
  else
  {
    // fp: derivative along the current descent direction p_
    // fp = df(x0+alpha*p)/dalpha at x
    fp = 0.0;
    for ( int i = 0; i < g.size(); i++ )
      fp += g[i] * p_[i];
    double new_alpha = linmin_.next_alpha(alpha_,f,fp);
    if ( debug_print )
      cout << "  CGOptimizer: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;

    if ( linmin_.fail() )
    {
      // line minimization failed
      if ( debug_print )
        cout << "CGOptimizer: line minimization failed" << endl;

      // restart from current point
      p_ = -g;

      x0_ = x;
      f0_ = f;

      g0norm2_ = norm2(g);
      if ( g0norm2_ < tol )
      {
        xp = x;
        return;
      }
      fp = -g0norm2_;
      fp0_ = fp;
      linmin_.reset();
      alpha_ = linmin_.next_alpha(alpha_,f,fp);
      if ( debug_print )
        cout << "  CGOptimizer: restart after fail: alpha=" << alpha_
           << " f=" << f << " fp=" << fp << endl;

      xp = x0_ + alpha_ * p_;
      first_step_ = false;
    }

    if ( linmin_.done() )
    {
      // wolfe1_ && wolfe2_ are true at alpha_
      if ( debug_print )
        cout << "  CGOptimizer: done with current descent direction" << endl;
      // define a new descent direction p_ using the Fletcher-Reeves formula
      assert(g0norm2_ > 0.0);
      double beta = norm2(g) / g0norm2_;
      if ( beta_max_ > 0.0 && beta > beta_max_ )
      {
        if ( debug_print )
          cout << "  CGOptimizer: beta exceeds beta_max " << endl;
        beta = min(beta, beta_max_);
      }
      if ( debug_print )
        cout << "  CGOptimizer: beta = " << beta << endl;
      p_ = beta * p_ - g;

      x0_ = x;
      f0_ = f;
      // recalculate f0, fp0
      // fp0 = d_e / d_alpha in direction pc_
      fp0_ = 0.0;
      for ( int i = 0; i < g.size(); i++ )
        fp0_ += g[i] * p_[i];
      g0norm2_ = norm2(g);
      fp = fp0_;
      assert(fp<0.0 && "CGOptimizer: p_ not a descent direction");

      // reset the line minimizer
      linmin_.reset();
      alpha_ = linmin_.next_alpha(alpha_,f,fp);
      if ( debug_print )
        cout << "  CGOptimizer: restart: alpha=" << alpha_
             << " f=" << f << " fp=" << fp << endl;
    }
    else
    {
      // not done with the current descent direction
      alpha_ = new_alpha;
    }
    xp = x0_ + alpha_ * p_;
  }
}
