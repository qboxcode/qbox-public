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
// CGOptimizer.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CGOPTIMIZER_H
#define CGOPTIMIZER_H

#include "LineMinimizer.h"
#include <valarray>

class CGOptimizer
{
  private:

  int n_;
  bool first_step_, debug_print;
  std::valarray<double> x0_, p_, gm_;
  double f0_, fp0_, g0norm2_, alpha_, beta_max_;
  LineMinimizer linmin_;
  double norm2(std::valarray<double>& v);

  public:

  CGOptimizer(int n): n_(n), first_step_(true), alpha_(0.0), beta_max_(0.0),
    debug_print(false)
  {
    x0_.resize(n);
    p_.resize(n);
    gm_.resize(n);
  }

  void reset(void) { first_step_ = true; }
  void set_sigma1(double s) { linmin_.set_sigma1(s); }
  void set_sigma2(double s) { linmin_.set_sigma2(s); }
  void set_alpha_start(double a ) { linmin_.set_alpha_start(a); }
  void set_alpha_max(double a ) { linmin_.set_alpha_max(a); }
  void set_beta_max(double b ) { beta_max_ = b; }
  void set_debug_print(void) { debug_print = true; linmin_.set_debug_print(); }

  int size(void) const { return n_; }
  double sigma1(void) const { return linmin_.sigma1(); }
  double sigma2(void) const { return linmin_.sigma2(); }
  double alpha(void) const { return alpha_; }
  double alpha_start(void) const { return linmin_.alpha_start(); }
  double beta_max(void) const { return beta_max_; }
  void compute_xp(const std::valarray<double>& x, const double f,
                  std::valarray<double>& g, std::valarray<double>& xp);
};
#endif
