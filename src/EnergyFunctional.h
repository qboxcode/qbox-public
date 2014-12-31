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
// EnergyFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ENERGYFUNCTIONAL_H
#define ENERGYFUNCTIONAL_H

#include <complex>
#include <vector>
#include <valarray>
#include <map>
#include <string>
#include "ChargeDensity.h"
#include "StructureFactor.h"
#include "ElectricEnthalpy.h"
#include "Timer.h"

class D3vector;
class Sample;
class Basis;
class AtomSet;
class Wavefunction;
class UnitCell;
class FourierTransform;
class XCOperator;
class NonLocalPotential;
class ConfinementPotential;

typedef std::map<std::string,Timer> TimerMap;

class EnergyFunctional
{
  private:

  Sample& s_;
  const ChargeDensity& cd_;
  Basis* vbasis_;
  FourierTransform *vft;
  std::vector<FourierTransform*> ft;
  StructureFactor sf;
  XCOperator* xco;
  ElectricEnthalpy* el_enth_;
  std::vector<std::vector<NonLocalPotential*> > nlp;    // nlp[ispin][ikp]
  std::vector<ConfinementPotential*> cfp; // cfp[ikp]

  std::vector<std::vector<double> > vps, dvps, rhops;
  std::vector<std::complex<double> > tmp_r, vion_local_g,
    dvion_local_g, vlocal_g, rhopst, rhogt, rhoelg, vtemp;
  std::vector<double> ftmp;

  std::vector<std::vector<double> > tau0, fion_esr;
  std::vector<std::vector<double> > fext;
  std::vector<double> zv_, rcps_;
  std::vector<int> na_;
  int namax_;
  int nsp_;
  double ekin_, econf_, eps_, enl_, ehart_,
         ecoul_, exc_, esr_, eself_, ets_, eexf_, exhf_, etotal_;
  double epv_, eefield_, enthalpy_;
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
  double exhf(void) const { return exhf_; }
  double esr(void) const { return esr_; }
  double eself(void) const { return eself_; }
  double ets(void) const { return ets_; }
  double eexf(void) const { return eexf_; }
  double eefield(void) const { return eefield_; }
  double epv(void) const { return epv_; }
  double enthalpy(void) const { return enthalpy_; }

  ElectricEnthalpy* el_enth() { return el_enth_; }

  const ConfinementPotential *confpot(int ikp) const { return cfp[ikp]; }

  void update_vhxc(void);

  void atoms_moved(void);
  void cell_moved(void);

  void print(std::ostream& os) const;

  EnergyFunctional(Sample& s, const ChargeDensity& cd);
  ~EnergyFunctional();
};
std::ostream& operator << ( std::ostream& os, const EnergyFunctional& e );
#endif
