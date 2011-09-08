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
// LineMinimizer.C
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
  assert(fp_low*(alpha_high-alpha_low)<0);
  // if fp_low * fp_high > 0, use bisection
  if ( fp_low * fp_high > 0 )
  {
    return 0.5 * ( alpha_low + alpha_high );
  }
  else
  {
    const double dalpha = alpha_high - alpha_low;
    return alpha_low - 0.5 * ( fp_low * dalpha * dalpha ) /
                             ( f_high - f_low - fp_low * dalpha );
  }
}
////////////////////////////////////////////////////////////////////////////////
double LineMinimizer::next_alpha(double alpha, double f, double fp)
{
  nstep_++;
  if ( nstep_ > nstep_max_ )
  {
#ifdef DEBUG
    cout << "LineMinimizer: fail, nstep_max" << endl;
#endif
    fail_ = true;
  }

  if ( done_ || fail_ )
    return alpha;

  if ( first_use )
  {
    first_use = false;
    f0 = f;
    fp0 = fp;
    alpha_m = 0;
    fm = f;
    fpm = fp;
    return alpha_start_;
  }

  // !first_use

  bool wolfe1 = f < f0 + sigma1_ * alpha * fp0;
  bool wolfe2 = fabs(fp) < sigma2_ * fabs(fp0);

  if ( !bracketing )
  {
    if ( !wolfe1 || f > fm )
    {
#ifdef DEBUG
      cout << "LineMinimizer: not bracketing, !wolfe1 || f > fm" << endl;
      if ( !wolfe1 ) cout << "LineMinimizer: !wolfe1" << endl;
      if ( f > fm ) cout << "LineMinimizer: f > fm" << endl;
#endif
      bracketing = true;
      alpha_low = alpha_m;
      f_low = fm;
      fp_low = fpm;
      alpha_high = alpha;
      f_high = f;
      fp_high = fp;
#ifdef DEBUG
      cout << "LineMinimizer: start bracketing" << endl;
      cout << "LineMinimizer: alpha=" << alpha << endl;
      cout << "LineMinimizer: alpha_low/alpha_high= "
           << alpha_low << " " << alpha_high << endl;
      cout << "LineMinimizer: fp_low=" << fp_low << endl;
#endif
      assert(fp_low*(alpha_high-alpha_low)<0);
      return interpolate();
    }

    // wolfe1 == true && f <= fm

    if ( wolfe2 )
    {
      // wolfe1 and wolfe2 both ok
      done_ = true;
      return alpha;
    }

    // wolfe1==true, wolfe2==false, f <= fm

    if ( fp >= 0 )
    {
#ifdef DEBUG
      cout << "LineMinimizer: not bracketing, fp>=0" << endl;
#endif
      bracketing = true;
      alpha_low = alpha;
      f_low = f;
      fp_low = fp;
      alpha_high = alpha_m;
      f_high = fm;
      fp_high = fpm;
      assert(fp_low*(alpha_high-alpha_low)<0);
      return interpolate();
    }

    // wolfe1==true, wolfe2==false, fp<0

    const double delta_alpha = 1.1 * (alpha - alpha_m);
    // increase alpha by at least 0.1*alpha_start
    const double new_alpha = alpha + std::max(delta_alpha,0.1*alpha_start_);
    alpha_m = alpha;
    fm = f;
    fpm = fp;
    return new_alpha;
  }
  else
  {
    // bracketing
    if ( !wolfe1 || f > f_low )
    {
#ifdef DEBUG
      cout << "LineMinimizer: bracketing, !wolfe1 || f > f_low" << endl;
      if ( !wolfe1 ) cout << "LineMinimizer: !wolfe1" << endl;
      if ( f > f_low ) cout << "LineMinimizer: f > f_low" << endl;
#endif
      alpha_high = alpha;
      f_high = f;
      fp_high = fp;
      assert(fp_low*(alpha_high-alpha_low)<0);
      return interpolate();
    }

    // wolfe1==true, f <= f_low

    if ( wolfe2 )
    {
      // both wolfe1 and wolfe2 ok
      done_ = true;
      return alpha;
    }

    // wolfe1==true, wolfe2==false, f<=f_low

    if ( fp*(alpha_high-alpha_low) >= 0 )
    {
#ifdef DEBUG
      cout << "LineMinimizer: bracketing, path1" << endl;
#endif
      alpha_high = alpha_low;
      f_high = f_low;
      fp_high = fp_low;
      alpha_low = alpha;
      f_low = f;
      fp_low = fp;
      assert(fp_low*(alpha_high-alpha_low)<0);
      return interpolate();
    }
    else
    {
      // there is an inflection point in [alpha_low,alpha_high]
      // this is likely caused by inaccuracy in f, fp or both
      // signal failure and exit
#ifdef DEBUG
      cout << "LineMinimizer: bracketing, path2 (fail)" << endl;
#endif
      fail_ = true;
      return alpha;
    }
  }
}
