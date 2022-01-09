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
// IonicStepper.cpp:
//
////////////////////////////////////////////////////////////////////////////////

#include "IonicStepper.h"

IonicStepper::IonicStepper (Sample& s) : s_(s), atoms_(s.atoms),
  constraints_(s.constraints), dt_(s.ctrl.dt)
{
  ndofs_ = 3 * atoms_.size() - constraints_.ndofs();
  // if there are more constraints than dofs, set ndofs_ to zero
  if ( ndofs_ < 0 ) ndofs_ = 0;
  nsp_ = atoms_.nsp();
  na_.resize(nsp_);
  r0_.resize(nsp_);
  rp_.resize(nsp_);
  v0_.resize(nsp_);
  pmass_.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ )
  {
    const int nais = atoms_.na(is);
    na_[is] = nais;
    r0_[is].resize(3*nais);
    rp_[is].resize(3*nais);
    v0_[is].resize(3*nais);
    pmass_[is] = atoms_.species_list[is]->mass() * 1822.89;
  }
  natoms_ = atoms_.size();
  // get positions and velocities from atoms_
  get_positions();
  get_velocities();
}

// center of mass position
D3vector IonicStepper::rcm(const std::vector<std::vector<double> >& r)
{
  D3vector mrsum;
  double msum = 0.0;
  for ( int is = 0; is < nsp_; is++ )
  {
    int nais = na_[is];
    double mass = pmass_[is];
    for ( int ia = 0; ia < nais; ia++ )
    {
      mrsum += mass * D3vector(r[is][3*ia],r[is][3*ia+1],r[is][3*ia+2]);
      msum += mass;
    }
  }
  if ( msum == 0.0 ) return D3vector(0,0,0);
  return mrsum / msum;
}

// reset center of mass position of r2 to be equal to that of r1
void IonicStepper::reset_rcm(const std::vector<std::vector<double> >& r1,
  std::vector<std::vector<double> >& r2)
{
  D3vector drcm = rcm(r2) - rcm(r1);
  for ( int is = 0; is < nsp_; is++ )
  {
    int nais = na_[is];
    for ( int ia = 0; ia < nais; ia++ )
    {
      r2[is][3*ia]   -= drcm.x;
      r2[is][3*ia+1] -= drcm.y;
      r2[is][3*ia+2] -= drcm.z;
    }
  }
}

// center of mass velocity
D3vector IonicStepper::vcm(const std::vector<std::vector<double> >& v)
{
  D3vector mvsum;
  double msum = 0.0;
  for ( int is = 0; is < nsp_; is++ )
  {
    int nais = na_[is];
    double mass = pmass_[is];
    for ( int ia = 0; ia < nais; ia++ )
    {
      mvsum += mass * D3vector(v[is][3*ia],v[is][3*ia+1],v[is][3*ia+2]);
      msum += mass;
    }
  }
  if ( msum == 0.0 ) return D3vector(0,0,0);
  return mvsum / msum;
}

// reset center of mass velocity to zero
void IonicStepper::reset_vcm(std::vector<std::vector<double> >& v)
{
  D3vector dvcm = vcm(v);
  for ( int is = 0; is < nsp_; is++ )
  {
    int nais = na_[is];
    for ( int ia = 0; ia < nais; ia++ )
    {
      v[is][3*ia]   -= dvcm.x;
      v[is][3*ia+1] -= dvcm.y;
      v[is][3*ia+2] -= dvcm.z;
    }
  }
}
