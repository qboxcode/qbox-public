////////////////////////////////////////////////////////////////////////////////
//
// IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: IonicStepper.h,v 1.1 2003-11-21 20:01:06 fgygi Exp $

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
  
  virtual void preprocess(const vector<vector<double> >&fion) {}
  virtual void postprocess(const vector<vector<double> >&fion) {}
  virtual void update(const vector<vector< double> >& fion) = 0;
  virtual double ekin(void) const = 0;
  virtual double temp(void) const = 0;
  
  virtual ~IonicStepper() {}
    
};
#endif
