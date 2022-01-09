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
// IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef IONICSTEPPER_H
#define IONICSTEPPER_H

#include "Sample.h"
#include "Species.h"
#include <vector>

class IonicStepper
{
  protected:

  Sample& s_;
  AtomSet& atoms_;
  ConstraintSet& constraints_;
  double dt_;
  int    nsp_;
  int    natoms_;
  // ndofs_ is the total number of degrees of freedom after
  // constraints are considered
  int                            ndofs_;
  std::vector<int>               na_;  // number of atoms per species na_[nsp_]
  std::vector<std::vector< double> >  r0_; // r0_[nsp_][3*na_]
  std::vector<std::vector< double> >  rp_; // rp_[nsp_][3*na_]
  std::vector<std::vector< double> >  rm_; // rm_[nsp_][3*na_]
  std::vector<std::vector< double> >  v0_; // v0_[nsp_][3*na_]
  std::vector<double>            pmass_;   // pmass_[nsp_]

  public:

  IonicStepper (Sample& s) : s_(s), atoms_(s.atoms),
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

  double r0(int is, int i) const { return r0_[is][i]; }
  double v0(int is, int i) const { return v0_[is][i]; }
  const std::vector<std::vector<double> >& r0(void) const { return r0_; }
  const std::vector<std::vector<double> >& v0(void) const { return v0_; }
  const std::vector<double>& pmass(void) const { return pmass_; }
  void get_positions(void) { atoms_.get_positions(r0_); }
  void get_velocities(void) { atoms_.get_velocities(v0_); }
  void set_positions(void) { atoms_.set_positions(r0_); }
  void set_velocities(void) { atoms_.set_velocities(v0_); }
  // center of mass position
  D3vector rcm(const std::vector<std::vector<double> >& r)
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
  void reset_rcm(const std::vector<std::vector<double> >& r1,
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
  D3vector vcm(const std::vector<std::vector<double> >& v)
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
  void reset_vcm(std::vector<std::vector<double> >& v)
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

  void setup_constraints(void)
  {
    constraints_.setup(atoms_);
  }
  virtual void compute_r(double e0,
    const std::vector<std::vector< double> >& f0) = 0;
  virtual void compute_v(double e0,
    const std::vector<std::vector< double> >& f0) = 0;
  virtual void reset(void) {}
  virtual double ekin(void) const { return 0.0; }
  virtual double ekin_stepper(void) const { return 0.0; }
  virtual double temp(void) const { return 0.0; }

  virtual ~IonicStepper() {}

};
#endif
