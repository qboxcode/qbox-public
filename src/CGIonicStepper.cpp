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
// CGIonicStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "CGIonicStepper.h"
#include "CGOptimizer.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CGIonicStepper::CGIonicStepper(Sample& s) : IonicStepper(s),
  cgopt_(CGOptimizer(3*natoms_))
{
  cgopt_.set_alpha_start(1.0);
  cgopt_.set_alpha_max(50.0);
  cgopt_.set_beta_max(10.0);
#ifdef DEBUG
  if ( MPIdata::onpe0() )
    cgopt_.set_debug_print();
#endif
}
////////////////////////////////////////////////////////////////////////////////
void CGIonicStepper::compute_r(double e0, const vector<vector<double> >& f0)
{
  // CG algorithm

  valarray<double> x(3*natoms_),xp(3*natoms_),g(3*natoms_);
  vector<vector<double> > gvec;
  gvec.resize(r0_.size());
  for ( int is = 0, i = 0; is < r0_.size(); is++ )
  {
    gvec[is].resize(r0_[is].size());
    for ( int j = 0; j < r0_[is].size(); j++ )
    {
      x[i] = r0_[is][j];
      gvec[is][j] = -f0[is][j];
      i++;
    }
  }

  // enforce compatibility of the gradient with constraints
  constraints_.enforce_v(r0_,gvec);

  // copy projected gradient to g
  for ( int is = 0, i = 0; is < r0_.size(); is++ )
    for ( int j = 0; j < r0_[is].size(); j++ )
      g[i++] = gvec[is][j];

  cgopt_.compute_xp(x,e0,g,xp);

  // check largest displacement
  // max_disp: largest acceptable displacement
  const double max_disp = 0.2;
  double largest_disp = 0.0;
  for ( int i = 0; i < xp.size(); i++ )
    largest_disp = max(largest_disp,fabs(xp[i]-x[i]));
  if ( largest_disp > max_disp )
  {
    if ( MPIdata::onpe0() )
      cout << "  CGIonicStepper: displacement exceeds limit, rescaling" << endl;
    // rescale displacement and reset the CG optimizer
    double fac = max_disp/largest_disp;
    xp = x + fac * (xp - x);
    cgopt_.set_alpha_start(fac*cgopt_.alpha_start());
    cgopt_.reset();
  }

  if ( MPIdata::onpe0() )
  {
    cout << "  CGIonicStepper: alpha = " << cgopt_.alpha() << endl;
  }

  for ( int is = 0, i = 0; is < r0_.size(); is++ )
    for ( int j = 0; j < r0_[is].size(); j++ )
      rp_[is][j] = xp[i++];

  if ( s_.ctrl.lock_cm )
    reset_rcm(r0_,rp_);

  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
  atoms_.reset_velocities();
}
