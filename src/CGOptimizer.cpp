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
#include <numeric>
#include "blas.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////////
void CGOptimizer::compute_xp(const valarray<double>& x, const double f,
                             valarray<double>& g, valarray<double>& xp)
{
  // Use the function value f and the gradient g at x to generate a new point xp
  // using the Fletcher-Reeves or Polak-Ribiere CG algorithm
  // return xp=x if the 2-norm of g is smaller than tol
  const double tol = 1.0e-18;
  const int one = 1;

  assert(x.size()==n_ && g.size()==n_ && xp.size()==n_);

  double fp;
  // define the descent direction
  if ( first_step_ )
  {
    p_ = -g;
    gm_ = g;

    x0_ = x;
    f0_ = f;

    g0norm2_ = ddot(&n_,&g[0],&one,&g[0],&one);
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
    fp = ddot(&n_,&g[0],&one,&p_[0],&one);
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

      g0norm2_ = ddot(&n_,&g[0],&one,&g[0],&one);
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
      return;
    }

    if ( linmin_.done() )
    {
      // wolfe1_ && wolfe2_ are true at alpha_
      if ( debug_print )
        cout << "  CGOptimizer: done with current descent direction" << endl;
      // define a new descent direction p_ using the Fletcher-Reeves formula
      assert(g0norm2_ > 0.0);

#if 0
      // Fletcher-Reeves
      double beta = ddot(&n_,&g[0],&one,&g[0],&one) / g0norm2_;
#else
      // Polak-Ribiere
      double beta = (ddot(&n_,&g[0],&one,&g[0],&one)-
                     ddot(&n_,&gm_[0],&one,&g[0],&one)) / g0norm2_;
#endif

      if ( beta_max_ > 0.0 && fabs(beta) > beta_max_ )
      {
        if ( debug_print )
          cout << "  CGOptimizer: |beta| exceeds beta_max " << endl;
        beta = (beta > 0.0) ? beta_max_ : -beta_max_;
      }
      if ( debug_print )
        cout << "  CGOptimizer: beta = " << beta << endl;
      p_ = beta * p_ - g;

      x0_ = x;
      f0_ = f;
      // recalculate f0, fp0
      // fp0 = d_e / d_alpha in direction pc_
      fp = ddot(&n_,&g[0],&one,&p_[0],&one);
      g0norm2_ = ddot(&n_,&g[0],&one,&g[0],&one);

      gm_ = g;
      fp0_ = fp;

      if ( fp > 0.0 )
      {
        // p_ is not a descent direction
        // restart from current point
        if ( debug_print )
          cout << "  CGOptimizer: p_ not a descent direction" << endl;
        p_ = -g;
        gm_ = g;

        x0_ = x;
        f0_ = f;

        g0norm2_ = ddot(&n_,&g[0],&one,&g[0],&one);
        fp = -g0norm2_;
        fp0_ = fp;
      }

      // set the starting alpha of the minimizer to be the current alpha_
      if ( alpha_ < linmin_.alpha_max() )
        linmin_.set_alpha_start(alpha_);
      else
        linmin_.set_alpha_start(0.5*linmin_.alpha_max());

      // reset the line minimizer
      linmin_.reset();
      alpha_ = linmin_.next_alpha(alpha_,f,fp);
      if ( debug_print )
        cout << "  CGOptimizer: restart: alpha=" << alpha_
             << " f=" << f << " fp=" << fp << endl;
    }
    xp = x0_ + alpha_ * p_;
  }
}
