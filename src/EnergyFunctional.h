////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.h,v 1.10 2004-02-04 19:55:17 fgygi Exp $

#ifndef ENERGYFUNCTIONAL_H
#define ENERGYFUNCTIONAL_H

#include <complex>
#include <vector>
#include <valarray>
#include <map>
#include <string>
#include "ChargeDensity.h"
#include "StructureFactor.h"
#include "Timer.h"
using namespace std;

class Sample;
class Basis;
class AtomSet;
class Wavefunction;
class UnitCell;
class FourierTransform;
class XCPotential;
class NonLocalPotential;

typedef map<string,Timer> TimerMap;

class EnergyFunctional
{
  private:
  
  const Sample& s_;
  ChargeDensity cd_;
  Basis* vbasis_;
  FourierTransform *vft;
  vector<FourierTransform*> ft;
  StructureFactor sf;
  XCPotential* xcp;
  NonLocalPotential* nlp;
  
  vector<vector<double> > vps, dvps, rhops;
  vector<complex<double> > tmp_r, vion_local_g, dvion_local_g, vlocal_g,
    rhopst, rhogt, rhoelg, vtemp;
  vector<double> ftmp;
  vector<vector<double> > v_r;
  
  vector<vector<double> > tau0, taum, fion_esr;
  vector<double> zv_, rcps_;
  vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, econf_, eps_, enl_, ehart_, 
         ecoul_, exc_, esr_, eself_, etotal_;
  valarray<double> fstress_, dfstress_;
  valarray<double> sigma_ekin,sigma_econf,sigma_eps,sigma_ehart,sigma_exc,
    sigma_enl, sigma_esr, sigma;

  void init(void);
  
  public:

  mutable TimerMap tmap;
  
  double energy(bool compute_hpsi, Wavefunction& dwf,
    bool compute_forces, vector<vector<double> >& fion,
    bool compute_stress, valarray<double>& sigma);
  
  double etotal(void) const { return etotal_; }
  double ekin(void) const { return ekin_; }
  double econf(void) const { return econf_; }
  double eps(void) const { return eps_; }
  double enl(void) const { return enl_; }
  double ehart(void) const { return ehart_; }
  double ecoul(void) const { return ecoul_; }
  double exc(void) const { return exc_; }
  double esr(void) const { return esr_; }
  double eself(void) const { return eself_; }
  
  void atoms_moved(void);
  void cell_moved(void);

  EnergyFunctional(const Sample& s);
  ~EnergyFunctional();
};
#endif
