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
// ANDIonicStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "ANDIonicStepper.h"
#include "AndersonMixer.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ANDIonicStepper::ANDIonicStepper(Sample& s) : IonicStepper(s),
  mixer_(3*natoms_,6,0)
{
  mixer_.restart();
}
////////////////////////////////////////////////////////////////////////////////
void ANDIonicStepper::compute_r(double e0, const vector<vector<double> >& f0)
{
  // Anderson mixer algorithm

  valarray<double> x(3*natoms_),xp(3*natoms_),g(3*natoms_);
  valarray<double> xbar(3*natoms_),fbar(3*natoms_);
  vector<vector<double> > gvec;
  gvec.resize(r0_.size());
  for ( int is = 0, i = 0; is < r0_.size(); is++ )
  {
    gvec[is].resize(r0_[is].size());
    for ( int j = 0; j < r0_[is].size(); j++ )
    {
      x[i] = r0_[is][j];
      gvec[is][j] = f0[is][j];
      i++;
    }
  }

  // enforce compatibility of the gradient with constraints
  constraints_.enforce_v(r0_,gvec);

  // copy projected gradient to g
  for ( int is = 0, i = 0; is < r0_.size(); is++ )
    for ( int j = 0; j < r0_[is].size(); j++ )
      g[i++] = gvec[is][j];

  if ( e0 > em_ ) mixer_.restart();
  em_ = e0;

  mixer_.update(&x[0],&g[0],&xbar[0],&fbar[0]);

  // check largest displacement
  // max_disp: largest acceptable displacement
  const double max_disp = 0.2;
  double largest_disp = 0.0;
  for ( int i = 0; i < xp.size(); i++ )
    largest_disp = max(largest_disp,fabs(xbar[i]-fbar[i]-x[i]));
  if ( largest_disp > max_disp )
  {
    if ( MPIdata::onpe0() )
      cout << "  ANDIonicStepper: displacement exceeds limit" << endl;
    // rescale displacement and reset the CG optimizer
    double fac = max_disp/largest_disp;
    xp = xbar + fac * fbar;
    mixer_.restart();
  }

  for ( int is = 0, i = 0; is < r0_.size(); is++ )
    for ( int j = 0; j < r0_[is].size(); j++ )
    {
      rp_[is][j] = xbar[i] + fbar[i];
      i++;
    }

  if ( s_.ctrl.lock_cm )
    reset_rcm(r0_,rp_);

  constraints_.enforce_r(r0_,rp_);
  rm_ = r0_;
  r0_ = rp_;
  atoms_.set_positions(r0_);
  atoms_.reset_velocities();
}
