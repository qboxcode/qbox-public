////////////////////////////////////////////////////////////////////////////////
//
// AtomSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomSet.h,v 1.10 2004-03-11 21:52:32 fgygi Exp $

#ifndef ATOMSET_H
#define ATOMSET_H

#include "Context.h"
#include "Atom.h"
#include "Species.h"
#include <list>
#include <map>
#include <string>
using namespace std;

class AtomSet
{
  private:
  
  const Context& ctxt_;
  
  int nel_;
  map<string,int> na_;  // na_[sp_name]: number of atoms of species sp_name
  map<string,int> isp_; // isp_[sp_name]: index of species sp_name
  vector<string> spname; // spname[is]: name of species is
  
  public:
  
  AtomSet(const Context& ctxt) : ctxt_(ctxt) {}

  vector<vector<Atom *> > atom_list; // atom_list[is][ia]
  vector<Species *> species_list;    // species_list[is]

  const Context& context(void) const { return ctxt_; }
  bool addAtom(Atom *a);
  bool delAtom(string name);
  bool addSpecies(Species *sp, string name);
  bool delSpecies(string name);
  bool reset(void); // remove all atoms and species
  Atom *findAtom(string name) const;
  Species *findSpecies(string name) const;
  void listAtoms(void) const;
  void listSpecies(void) const;
  int na(string spname) const;  // number of atoms of species spname
  int na(int is) const;         // number of atoms of species is
  int isp(string spname) const; // index of species spname
  int nel(void) const { return nel_; };
  int nsp(void) const { return species_list.size(); }
  void get_positions(vector<vector<double> >& tau) const;
  void set_positions(const vector<vector<double> >& tau);
  void get_velocities(vector<vector<double> >& vel) const;
  void set_velocities(const vector<vector<double> >& vel);
  void reset_velocities(void);
  int size(void);
};
ostream& operator << ( ostream &os, AtomSet &as );
#endif
