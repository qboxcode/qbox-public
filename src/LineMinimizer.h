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
// LineMinimizer.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LINEMINIMIZER_H
#define LINEMINIMIZER_H

class LineMinimizer
{
  private:

  double f0,fp0,fm,fpm,f_low,fp_low,f_high,fp_high,
         alpha_m,alpha_low,alpha_high;
  bool first_use, done_, fail_, bracketing;
  double alpha_start_, sigma1_, sigma2_;
  int nstep_, nstep_max_;

  double interpolate(void);

  public:

  LineMinimizer(void) : sigma1_(0.1), sigma2_(0.5), alpha_start_(1.0), 
   first_use(true), done_(false), fail_(false), bracketing(false), nstep_(0),
   nstep_max_(10) {}
  void reset(void) { first_use = true; done_ = false; fail_ = false;
    bracketing = false; nstep_ = 0; }
  double sigma1(void) const { return sigma1_; }
  double sigma2(void) const { return sigma2_; }
  double alpha_start(void) const { return alpha_start_; }
  bool done(void) const { return done_; }
  bool fail(void) const { return fail_; }

  void set_sigma1(double s) { sigma1_ = s; }
  void set_sigma2(double s) { sigma2_ = s; }
  void set_alpha_start(double a) { alpha_start_ = a; }

  double next_alpha(double alpha, double f, double fp);
};
#endif
