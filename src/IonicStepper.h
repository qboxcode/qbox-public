////////////////////////////////////////////////////////////////////////////////
//
// IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: IonicStepper.h,v 1.2 2003-11-27 01:18:41 fgygi Exp $

#ifndef IONICSTEPPER_H
#define IONICSTEPPER_H

#include "Sample.h"
#include <vector>
using namespace std;

class IonicStepper
{
  protected:
  
  AtomSet& atoms_;
  double                    dt_;
  int                       nsp_;
  int                       ndofs_;
  vector<int>               na_;      // number of atoms per species na_[nsp_]
  vector<vector< double> >  tau0_;    // tau0_[nsp_][3*na_]
  vector<vector< double> >  vel_;     // vel_[nsp_][3*na_]
  vector<double>            pmass_;   // pmass_[nsp_]
  double                    ekin_;    // kinetic energy

  public:
  
  IonicStepper (Sample& s) : atoms_(s.atoms), dt_(s.ctrl.dt)
  {
    ndofs_ = 3 * s.atoms.size();
    nsp_ = atoms_.nsp();
    na_.resize(nsp_);
    tau0_.resize(nsp_);
    vel_.resize(nsp_);
    pmass_.resize(nsp_);
    for ( int is = 0; is < nsp_; is++ )
    {
      const int nais = atoms_.na(is);
      na_[is] = nais;
      tau0_[is].resize(3*nais);
      vel_[is].resize(3*nais);
      pmass_[is] = atoms_.species_list[is]->mass() * 1822.89;
    }
    atoms_.get_positions(tau0_);
    atoms_.get_velocities(vel_);
  }
  double tau0(int is, int i) const { return tau0_[is][i]; }
  double vel(int is, int i) const { return vel_[is][i]; }
  const vector<vector<double> >& tau0(void) const { return tau0_; }
  
  virtual void preprocess(const vector<vector<double> >&fion) {}
  virtual void postprocess(const vector<vector<double> >&fion) {}
  virtual void update(const vector<vector< double> >& fion) = 0;
  double ekin(void) const { return ekin_; }
  double temp(void) const
  {
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    if ( ndofs_ > 0.0 )
      return 2.0 * ( ekin_ / boltz ) / ndofs_;
    else
      return 0.0;
  }
  
  virtual ~IonicStepper() {}
    
};
#endif
