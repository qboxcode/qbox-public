////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.h,v 1.9 2003-06-11 22:10:11 fgygi Exp $

#ifndef ENERGYFUNCTIONAL_H
#define ENERGYFUNCTIONAL_H

#include <complex>
#include <vector>
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
  vector<complex<double> > tmp_r, vion_local_g, vlocal_g, rhopst,
    rhogt, rhoelg, vtemp;
  vector<double> ftmp;
  vector<vector<double> > v_r;
  
  vector<vector<double> > tau0, taum, fion_esr, fion;
  vector<double> zv_, rcps_;
  vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, eps_, enl_, ehart_, ecoul_, exc_, esr_, eself_, etotal_;

  void init(void);
  
  public:

  mutable TimerMap tmap;
  
  double energy(bool compute_hpsi, Wavefunction& dwf,
    bool compute_forces, vector<vector<double> >& fion,
    bool compute_stress, UnitCell& dcell);
  
  double etotal(void) const { return etotal_; }
  double ekin(void) const { return ekin_; }
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
