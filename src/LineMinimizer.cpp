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
// LineMinimizer.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "LineMinimizer.h"
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;
////////////////////////////////////////////////////////////////////////////////
double LineMinimizer::interpolate(void)
{
  const double dalpha = alpha_high - alpha_low;
  // use psi'(alpha_low), psi(alpha_low_), psi(alpha_high)

  double new_alpha;
  if ( use_psi )
  {
    // use psi
    double psip_low = psip(fp_low);
    double psi_low = psi(alpha_low,f_low);
    double psi_high = psi(alpha_high,f_high);
    new_alpha = alpha_low - 0.5 * ( psip_low * dalpha * dalpha ) /
                ( psi_high - psi_low - psip_low * dalpha );
  }
  else
  {
    // use f
    // quadratic interpolation using f_low, fp_low, f_high
    // new_alpha = alpha_low - 0.5 * ( fp_low * dalpha * dalpha ) /
    //                  ( f_high - f_low - fp_low * dalpha );
    if ( fp_low*fp_high < 0 )
    {
      // secant
      new_alpha = alpha_low - fp_low *
        ( (alpha_high-alpha_low)/(fp_high-fp_low) );
    }
    else
    {
      // midpoint
      new_alpha = 0.5 * (alpha_low+alpha_high);
    }
  }

  if ( debug_print )
  {
    cout << "LineMinimizer: interpolate: [alpha_low,alpha_high] = ["
         << alpha_low << "," << alpha_high << "]" << endl;
    cout << "LineMinimizer: interpolate: f_low, f_high: "
         << f_low << " " << f_high << endl;
    cout << "LineMinimizer: interpolate: fp_low, fp_high: "
         << fp_low << " " << fp_high << endl;
    cout << "LineMinimizer: interpolate: new_alpha: " << new_alpha << endl;
  }
  if ( alpha_low < alpha_high )
  {
    assert ( new_alpha >= alpha_low && new_alpha <= alpha_high );
  }
  else
  {
    assert ( new_alpha <= alpha_low && new_alpha >= alpha_high );
  }
  return new_alpha;
}
////////////////////////////////////////////////////////////////////////////////
double LineMinimizer::next_alpha(double alpha, double f, double fp)
{
  if ( done_ || fail_ )
    return alpha;

  if ( first_use )
  {
    first_use = false;
    f0 = f;
    fp0 = fp;
    alpha_low = 0;
    f_low = f0;
    fp_low = fp0;
    alpha_high = alpha_max_;
    assert(alpha_max_ > alpha_start_);
    if ( debug_print )
      cout << "LineMinimizer: first use: f0, fp0: " << f0 << " " << fp0 << endl;
    return alpha_start_;
  }
  if ( debug_print )
    cout << "LineMinimizer: next_alpha(" << alpha << ","
         << f << "," << fp << ")" << endl;

  bool wolfe1 = f < f0 + sigma1_ * alpha * fp0;
  bool wolfe2 = fabs(fp) < sigma2_ * fabs(fp0);

  if ( debug_print )
  {
    cout << "LineMinimizer: wolfe1: f = " << f << endl;
    cout << "LineMinimizer: wolfe1: f0 + sigma1_ * alpha * fp0 = "
         << f0 + sigma1_ * alpha * fp0 << endl;
    cout << "LineMinimizer: wolfe1/wolfe2: " << wolfe1 << "/" << wolfe2 << endl;
  }

  // check if alpha satisfies both wolfe1 and wolfe2 and return
  if ( wolfe1 && wolfe2 )
  {
    done_ = true;
    return alpha;
  }

  if ( !bracketing )
  {
    // we have not entered the bracketing phase yet
    // Enter bracketing mode if condition U1 holds: psi(alpha) > psi(alpha_low)
    // Note: U1 is equivalent to wolfe1(alpha) == false
    if ( psi(alpha,f) > psi(alpha_low,f_low) )
    {
      // wolfe1(alpha) == false
      // we can enter the bracketing phase
      // enter the bracketing phase with alpha_low = 0, alpha_high = alpha

      if ( debug_print )
      {
        cout << "LineMinimizer: entering bracketing: wolfe1==false" << endl;
        cout << "LineMinimizer: psi(alpha), psip(alpha):"
             << psi(alpha,f) << " " << psip(fp) << endl;
        cout << "LineMinimizer: psi(alpha_low, psi(alpha_high):"
             << psi(alpha_low,f_low) << " " << psi(alpha_high,f_high) << endl;
        cout << "LineMinimizer: psip(alpha_low, psip(alpha_high):"
           << psip(fp_low) << " " << psip(fp_high) << endl;
      }

      bracketing = true;
      use_psi = true;
      alpha_high = alpha;
      f_high = f;
      fp_high = fp;
      assert(psi(alpha_low,f_low)<=0);
      assert(psip(fp_low)*(alpha_high-alpha_low)<=0);
      assert(psi(alpha_low,f_low)<=psi(alpha_high,f_high));
      return interpolate();
    }

    // check if wolfe1 is ok and f'(alpha)(alpha-alpha_low)  > 0
    // Need to test only the second inequality at this point
    if ( fp*(alpha-alpha_low) > 0 )
    {
      // enter bracketing mode with alpha_high = alpha_low, alpha_low = alpha
      if ( debug_print )
        cout << "LineMinimizer: entering bracketing: case U3" << endl;

      bracketing = true;
      use_psi = false;
      alpha_high = alpha_low;
      f_high = f_low;
      fp_high = fp_low;
      alpha_low = alpha;
      f_low = f;
      fp_low = fp;
      assert(psi(alpha_low,f_low)<=0);
      assert(psip(fp_low)*(alpha_high-alpha_low)<=0);
      assert(psi(alpha_low,f_low)<=psi(alpha_high,f_high));
      return interpolate();
    }

    // Condition U2 holds
    // increase alpha

    if ( debug_print )
      cout << "LineMinimizer: U2, increase alpha" << endl;

    double new_alpha = std::min(alpha+delta_*(alpha-alpha_low), alpha_max_);
    if ( new_alpha == alpha_max_ )
      done_ = true;
    return new_alpha;
  }
  else
  {
    // we are already in bracketing mode
    nstep_++;
    if ( nstep_max_ > 0 && nstep_ > nstep_max_ )
    {
      if ( debug_print )
        cout << "LineMinimizer: fail, nstep_max" << endl;

      fail_ = true;
      return alpha;
    }

    if ( debug_print )
    {
      cout << "LineMinimizer: bracketing mode: [alpha_low,alpha_high] = ["
           << alpha_low << "," << alpha_high << "]" << endl;
      cout << "LineMinimizer: bracketing mode: f_low, f_high: "
           << f_low << " " << f_high << endl;
      cout << "LineMinimizer: bracketing mode: fp_low, fp_high: "
           << fp_low << " " << fp_high << endl;
      cout << "LineMinimizer: bracketing mode: psi(alpha), psip(alpha):"
           << psi(alpha,f) << " " << psip(fp) << endl;
      cout << "LineMinimizer: bracketing mode: psi(alpha_low, psi(alpha_high):"
           << psi(alpha_low,f_low) << " " << psi(alpha_high,f_high) << endl;
      cout << "LineMinimizer: bracketing mode: psip(alpha_low, psip(alpha_high):"
           << psip(fp_low) << " " << psip(fp_high) << endl;
    }

    // check U1: psi(alpha) > psi(alpha_low)
    if ( psi(alpha,f) > psi(alpha_low,f_low) )
    {
      if ( debug_print )
        cout << "LineMinimizer: bracketing, U1" << endl;
      alpha_high = alpha;
      f_high = f;
      fp_high = fp;
      assert(psi(alpha_low,f_low)<=0);
      assert(psip(fp_low)*(alpha_high-alpha_low)<=0);
      assert(psi(alpha_low,f_low)<=psi(alpha_high,f_high));
      return interpolate();
    }
    else
    {
      // at this point psi(alpha) <= psi(alpha_low)
      // test condition U2: psi'(alpha)*(alpha_low-alpha) > 0
      if ( psip(fp)*(alpha_low-alpha) > 0 )
      {
        if ( debug_print )
          cout << "LineMinimizer: bracketing, U2" << endl;
        alpha_low = alpha;
        f_low = f;
        fp_low = fp;
        assert(psi(alpha_low,f_low)<=0);
        assert(psip(fp_low)*(alpha_high-alpha_low)<=0);
        assert(psi(alpha_low,f_low)<=psi(alpha_high,f_high));
        return interpolate();
      }
      else
      {
        if ( debug_print )
          cout << "LineMinimizer: bracketing, U3" << endl;
        alpha_high = alpha_low;
        f_high = f_low;
        fp_high = fp_low;
        alpha_low = alpha;
        f_low = f;
        fp_low = fp;
        assert(psi(alpha_low,f_low)<=0);
        assert(psip(fp_low)*(alpha_high-alpha_low)<=0);
        assert(psi(alpha_low,f_low)<=psi(alpha_high,f_high));
        use_psi = false;
        return interpolate();
      }
    }
  }
}
