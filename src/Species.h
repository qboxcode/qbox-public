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
// Species.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SPECIES_H
#define SPECIES_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>

class Species
{
  private:

  struct ProjectorData
  {
    int l, m, n;
  };

  enum PP_type
  {
    // norm conserving pseudopotential
    NCPP,
    // norm conserving semilocal pseudopotential
    SLPP
  };

  PP_type type_; // identify type of the PP

  int nlm_;             // number of non-local projectors (including m)
  int nop_;             // number of non-local projectors (excluding m)
  int ndft_;

  std::vector<std::vector<double> > vps_spl_, vps_spl2_, phi_spl_, phi_spl2_;
  std::vector<double>               gspl_, vlocg_spl_, vlocg_spl2_;
  std::vector<std::vector<double> > vnlg_spl_, vnlg_spl2_;
  std::vector<double> wsg_;  // wsg_[l] Kleinman-Bylander weight
                             // 1/<phi|delta_V|phi>
  std::vector<double>               nlcc_spl_, nlcc_spl2_;
  std::vector<double>               nlccg_spl_, nlccg_spl2_;

  std::vector<double> rps_spl_;  // radial linear mesh (same for all l)

  std::string name_;         // name used in current application
  std::string uri_;          // uri of the resource defining the pseudopotential

  std::string symbol_;
  int atomic_number_;
  double mass_;        // mass in a.m.u (Carbon = 12.0)

  std::string description_; // description of the pseudopotential
  int zval_;           // valence charge
  double zcore_;       // core charge
  double ztot_;        // total charge
  int lmax_;           // largest angular momentum
  int llocal_;         // angular momentum taken as local
  int nquad_;          // number of semi-local quadrature points
  double rquad_;       // end of semi-local quadrature interval
  double deltar_;      // mesh spacing for potentials and wavefunctions
  std::vector<std::vector<double> > vps_;  // potentials for each l (input)
  std::vector<std::vector<double> > phi_;  // atomic wf for each l (input)
  double rcps_;        // cutoff radius of gaussian pseudocharge
  std::vector<double> nlcc_; // non linear core correction

  // map index of projector -> angular momentum
  std::vector<int> lmap_;
  // local potential in radial representation
  std::vector<double> vlocal_;
  // projector in radial representation, indices proj_[l][n][r]
  std::vector<std::vector<std::vector<double> > > proj_;
  // matrix D in block diagonal storage, indices d_[l,n,m]
  std::vector<std::vector<std::vector<double> > > d_;

  // initialize a norm conserving PP
  bool initialize_ncpp();
  // initialize a semilocal potential
  bool initialize_slpp();
  // non linear core correction
  void initialize_nlcc();

  // helper function that extracts l, m and n from projector index
  ProjectorData get_proj_data(int ipr);

  public:

  Species(std::string name);

  const std::string& name(void) const { return name_; }
  const std::string& symbol(void) const { return symbol_; }
  const std::string& description(void) const { return description_; }
  const std::string& uri(void) const { return uri_; }
  int atomic_number(void) const { return atomic_number_; }
  int zval(void) const { return zval_; }
  double zcore(void) const { return zcore_; }
  double ztot(void) const { return ztot_; }
  double mass(void) const { return mass_; }
  int lmax(void) const { return lmax_; }
  int llocal(void) const { return llocal_; }
  int nquad(void) const { return nquad_; }
  double rquad(void) const { return rquad_; }
  double deltar(void) const { return deltar_; }
  double rcps(void) const { return rcps_; }
  bool has_nlcc(void) const { return nlcc_.size() > 0; }
  bool has_dmatrix(void) const { return d_.size() > 0; }

  // number of non-local projectors sum_(l!=llocal) (2*l+1)
  int nlm(void) { return nlm_; }
  // number of non-local projectors w/o m degeneracy
  int nop(void) { return nop_; }
  // angular momentum of projector with index iop
  int l(int iop) { return lmap_[iop]; }
  // extract D matrix
  double dmatrix(int ipr, int jpr);

  bool non_local(void) { return nop_ > 0; };
  double eself(void)
  { return ztot_ * ztot_ / ( sqrt ( 2.0 * M_PI ) * rcps_ ); };

  void phi(int l, double r, double &val);              // phi(l,r) in r space
  void vpsr(int l, double r, double &v);               // Vps(l,r) in r space
  void dvpsr(int l, double r, double &v, double &dv);  // Vps and dVps/dr
  void vlocg(double q, double &v);                     // Vloc(g)
  void dvlocg(double q, double &v, double &dv);        // Vloc(g) and dVloc/dg
  void vnlg(int iop, double q, double &v);             // Vnl(iop,g)
  void dvnlg(int iop, double q, double &v, double &dv);// Vnl(iop,g) and dVnl/dg
  double rhopsg(double q);        // pseudocharge in g space
  void rhocoreg(double q, double &rho);      // core correction in g space
  // core correction and its derivative in g space
  void drhocoreg(double q, double &rho, double &drho);

  double wsg(int iop) { return wsg_[iop]; };
  double rcut_loc(double epsilon); // radius beyond which potential is local

  const std::vector<std::vector<double> >& vps(void) const { return vps_spl_; }
  const std::vector<std::vector<double> >& phi(void) const { return phi_spl_; }

  bool initialize(double rcps);
  void info(std::ostream& os);
  void print(std::ostream& os, bool expanded_form);
  void print_nlcc(std::ostream& os);
  void print_radial_function(std::ostream& os, const std::vector<double>& rad_func);

  friend class SpeciesReader;
  friend class SpeciesHandler;

};
std::ostream& operator << ( std::ostream &os, Species &a );
class SpeciesInitException
{
  public:
  std::string msg;
  SpeciesInitException(std::string s) : msg(s) {}
};
#endif
