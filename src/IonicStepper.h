////////////////////////////////////////////////////////////////////////////////
//
// IonicStepper.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: IonicStepper.h,v 1.6 2004-12-17 23:31:00 fgygi Exp $

#ifndef IONICSTEPPER_H
#define IONICSTEPPER_H

#include "Sample.h"
#include <vector>
using namespace std;

class IonicStepper
{
  protected:
  
  Sample& s_;
  AtomSet& atoms_;
  double                    dt_;
  int                       nsp_;
  int                       ndofs_;
  vector<int>               na_;      // number of atoms per species na_[nsp_]
  vector<vector< double> >  r0_;      // r0_[nsp_][3*na_]
  vector<vector< double> >  rp_;      // rp_[nsp_][3*na_]
  vector<vector< double> >  v0_;      // v0_[nsp_][3*na_]
  vector<double>            pmass_;   // pmass_[nsp_]

  public:
  
  IonicStepper (Sample& s) : s_(s), atoms_(s.atoms), dt_(s.ctrl.dt)
  {
    ndofs_ = 3 * s.atoms.size();
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
    atoms_.get_positions(r0_);
  }
  
  double r0(int is, int i) const { return r0_[is][i]; }
  double rp(int is, int i) const { return rp_[is][i]; }
  double v0(int is, int i) const { return v0_[is][i]; }
  const vector<vector<double> >& r0(void) const { return r0_; }
  const vector<vector<double> >& rp(void) const { return rp_; }
  const vector<vector<double> >& v0(void) const { return v0_; }
  
  virtual void compute_rp(const vector<vector< double> >& f0) = 0;
  virtual void compute_v0(const vector<vector< double> >& f0) = 0;
  virtual void update_r(void) = 0;
  virtual void update_v(void) = 0;
  virtual void reset(void) {}
  virtual double ekin(void) const { return 0.0; }
  virtual double temp(void) const { return 0.0; }
  
  virtual ~IonicStepper() {}
    
};
#endif
