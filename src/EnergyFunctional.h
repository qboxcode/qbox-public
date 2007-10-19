////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.h,v 1.17 2007-10-19 17:10:58 fgygi Exp $

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

class Sample;
class Basis;
class AtomSet;
class Wavefunction;
class UnitCell;
class FourierTransform;
class XCPotential;
class NonLocalPotential;
class ConfinementPotential;

typedef std::map<std::string,Timer> TimerMap;

class EnergyFunctional
{
  private:

  const Sample& s_;
  const ChargeDensity& cd_;
  Basis* vbasis_;
  FourierTransform *vft;
  std::vector<FourierTransform*> ft;
  StructureFactor sf;
  XCPotential* xcp;
  std::vector<NonLocalPotential*> nlp;    // nlp[ikp]
  std::vector<ConfinementPotential*> cfp; // cfp[ikp]

  std::vector<std::vector<double> > vps, dvps, rhops;
  std::vector<std::complex<double> > tmp_r, vion_local_g,
    dvion_local_g, vlocal_g, rhopst, rhogt, rhoelg, vtemp;
  std::vector<double> ftmp;

  std::vector<std::vector<double> > tau0, fion_esr;
  std::vector<double> zv_, rcps_;
  std::vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, econf_, eps_, enl_, ehart_,
         ecoul_, exc_, esr_, eself_, ets_, etotal_;
  std::valarray<double> sigma_ekin,sigma_econf,sigma_eps,sigma_ehart,sigma_exc,
    sigma_enl, sigma_esr, sigma;

  public:

  std::vector<std::vector<double> > v_r;
  mutable TimerMap tmap;

  double energy(bool compute_hpsi, Wavefunction& dwf,
    bool compute_forces, std::vector<std::vector<double> >& fion,
    bool compute_stress, std::valarray<double>& sigma);

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
  double ets(void) const { return ets_; }

  const ConfinementPotential *confpot(int ikp) const { return cfp[ikp]; }

  void update_vhxc(void);

  void atoms_moved(void);
  void cell_moved(void);

  void print(std::ostream& os) const;

  EnergyFunctional(const Sample& s, const ChargeDensity& cd);
  ~EnergyFunctional();
};
std::ostream& operator << ( std::ostream& os, const EnergyFunctional& e );
#endif
