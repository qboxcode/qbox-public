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
// AtomSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomSet.h,v 1.24 2009-10-06 06:29:05 fgygi Exp $

#ifndef ATOMSET_H
#define ATOMSET_H

#include "Context.h"
#include "Atom.h"
#include "UnitCell.h"
#include <vector>
#include <string>
#include <list>
#include <map>
#include <string>

class Species;

class AtomSet
{
  private:

  const Context& ctxt_;

  int nel_;
  std::map<std::string,int> na_;  // na_[sp_name]: num. of at. of spec. sp_name
  std::map<std::string,int> isp_; // isp_[sp_name]: index of species sp_name
  std::map<std::string,int> is_; // is_[atom_name]: is index of atom atom_name
  std::map<std::string,int> ia_; // ia_[atom_name]: ia index of atom atom_name
  std::vector<std::string> spname; // spname[is]: name of species is

  UnitCell cell_;

  public:

  AtomSet(const Context& ctxt) : ctxt_(ctxt), nel_(0) {}
  ~AtomSet(void);

  std::vector<std::vector<Atom *> > atom_list; // atom_list[is][ia]
  std::vector<Species *> species_list;    // species_list[is]

  const Context& context(void) const { return ctxt_; }
  bool addAtom(Atom *a);
  bool delAtom(std::string name);
  bool addSpecies(Species *sp, std::string name);
  bool delSpecies(std::string name);
  bool reset(void); // remove all atoms and species
  Atom *findAtom(std::string name) const;
  Species *findSpecies(std::string name) const;
  void listAtoms(void) const;
  void listSpecies(void) const;
  int na(const std::string& spname) const;  // number of atoms of species spname
  int na(int is) const;         // number of atoms of species is
  int isp(const std::string& spname) const; // index of species spname
  int is(const std::string& atom_name) const;
  int ia(const std::string& atom_name) const;
  int nel(void) const { return nel_; };
  int nsp(void) const { return species_list.size(); }
  void get_positions(std::vector<std::vector<double> >& tau) const;
  void set_positions(const std::vector<std::vector<double> >& tau);
  void get_velocities(std::vector<std::vector<double> >& vel) const;
  void set_velocities(const std::vector<std::vector<double> >& vel);
  const UnitCell& cell(void) const { return cell_; }
  void set_cell(const UnitCell& cell) { cell_ = cell; }
  void set_cell(const D3vector& a, const D3vector& b, const D3vector& c)
  {
    cell_.set(a,b,c);
  }
  void sync(void);
  void reset_velocities(void);
  D3vector vcm(void) const;
  D3vector dipole(void) const;
  void reset_vcm(void);
  void fold_in_ws(void);
  int size(void) const;
 };
std::ostream& operator << ( std::ostream &os, const AtomSet &as );
#endif
