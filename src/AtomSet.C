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
// AtomSet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomSet.C,v 1.29 2010-04-16 21:40:50 fgygi Exp $

#include "AtomSet.h"
#include "Species.h"
#include "NameOf.h"
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AtomSet::~AtomSet(void)
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      delete atom_list[is][ia];
    }
    delete species_list[is];
  }
}
////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addSpecies(Species* sp, string name)
{
  const double rcps = 1.5;
  sp->initialize(rcps);

  Species *s = findSpecies(name);
  if ( s != 0 )
  {
    // species is already defined: substitute with new definition
    if ( ctxt_.onpe0() )
    {
      cout << " AtomSet::addSpecies: species " << name
           << " is already defined" << endl;
      cout << " AtomSet::addSpecies: redefining species" << endl;
    }
    // Check if s and sp are compatible: zval must be equal
    if ( s->zval() != sp->zval() )
    {
      cout << " AtomSet::addSpecies: species do not have the same"
           << " number of valence electrons. Cannot redefine species" << endl;
      return false;
    }

    int is = isp_[s->name()];
    species_list[is] = sp;
    delete s;
  }
  else
  {
    // create new entry in species list
    species_list.push_back(sp);
    spname.push_back(name);
    isp_[name] = species_list.size()-1;
    na_[name] = 0;
    atom_list.resize(atom_list.size()+1);
  }

  if ( ctxt_.onpe0() )
  {
    cout << endl << " species " << sp->name() << ":" << endl;
    sp->info(cout);
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addAtom(Atom *a)
{

  // check atom_list for name
  if ( findAtom(a->name()) )
  {
    // this name is already in the atom list, reject atom definition
    if ( ctxt_.onpe0() )
      cout << " AtomSet:addAtom: atom " << a->name()
           << " is already defined" << endl;
    return false;
  }

  // check if species is defined
  string spname = a->species();
  Species *s = findSpecies(spname);
  if ( !s )
  {
    // species not found, cannot define atom
    if ( ctxt_.onpe0() )
      cout << " AtomSet:addAtom: species " << spname
           << " is undefined" << endl;
    return false;
  }

  // add an atom to the atom_list
  int is = isp(spname);
  assert ( is >= 0 );
  atom_list[is].push_back(a);
  ia_[a->name()] = atom_list[is].size()-1;
  is_[a->name()] = is;

  // update count of atoms of species spname
  // increment count
  na_[spname]++;

  // update total number of electrons
  nel_ = 0;
  for ( int is = 0; is < species_list.size(); is++ )
  {
    nel_ += species_list[is]->zval() * atom_list[is].size();
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::delAtom(string name)
{
  vector<vector<Atom*> >::iterator psa = atom_list.begin();
  vector<Species*>::iterator ps = species_list.begin();
  for ( int is = 0; is < species_list.size(); is++ )
  {
    vector<Atom*>::iterator pa = atom_list[is].begin();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( atom_list[is][ia]->name() == name )
      {
        string spname = atom_list[is][ia]->species();
        na_[spname]--;
        delete atom_list[is][ia];

        // remove map entries ia_[name] and is_[name]
        map<string,int>::iterator i = ia_.find(name);
        ia_.erase(i);
        i = is_.find(name);
        is_.erase(i);

        atom_list[is].erase(pa);
        nel_ -= species_list[is]->zval();

        return true;
      }
      pa++;
    }
    psa++;
    ps++;
  }

  // this name was not found in the atom list
  if ( ctxt_.onpe0() )
    cout << " AtomSet:delAtom: no such atom: " << name << endl;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
Atom *AtomSet::findAtom(string name) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( atom_list[is][ia]->name() == name )
        return atom_list[is][ia];
    }
  }
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
Species *AtomSet::findSpecies(string name) const
{
  int is = isp(name);
  if ( is >= 0 )
    return species_list[is];
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::isp(const string& name) const
{
  map<string,int>::const_iterator i = isp_.find(name);
  if ( i != isp_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::is(const string& atom_name) const
{
  map<string,int>::const_iterator i = is_.find(atom_name);
  if ( i != is_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::ia(const string& atom_name) const
{
  map<string,int>::const_iterator i = ia_.find(atom_name);
  if ( i != ia_.end() )
    return (*i).second;
  else
    return -1;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::listAtoms(void) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      if ( ctxt_.onpe0() )
        cout << *atom_list[is][ia];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::listSpecies(void) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    if ( ctxt_.onpe0() )
    {
      cout << endl << " species " << spname[is] << ":" << endl;
      species_list[is]->info(cout);
      cout << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::na(const string& spname) const
{
  map<string,int>::const_iterator i = na_.find(spname);
  if ( i != na_.end() )
    return (*i).second;
  else
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::na(int is) const
{
  assert( is >= 0 && is < atom_list.size() );
  return atom_list[is].size();
}

////////////////////////////////////////////////////////////////////////////////
int AtomSet::size(void) const
{
  int n = 0;
  for ( int is = 0; is < atom_list.size(); is++ )
    n += atom_list[is].size();
  return n;
}

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::reset(void)
{
  // delete all atoms and species

  for ( int is = 0; is < species_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      delete atom_list[is][ia];
    }
    atom_list[is].resize(0);

    delete species_list[is];
  }
  atom_list.resize(0);
  species_list.resize(0);
  spname.resize(0);
  for ( map<string,int>::iterator i = na_.begin();
        i != na_.end();
        i++ )
        na_.erase(i);

  for ( map<string,int>::iterator i = isp_.begin();
        i != isp_.end();
        i++ )
        isp_.erase(i);

  for ( map<string,int>::iterator i = is_.begin();
        i != is_.end();
        i++ )
        is_.erase(i);

  for ( map<string,int>::iterator i = ia_.begin();
        i != ia_.end();
        i++ )
        ia_.erase(i);
  nel_ = 0;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_positions(vector<vector<double> >& tau) const
{
  if (tau.size() != atom_list.size())
    tau.resize(atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    if (tau[is].size() != 3*atom_list[is].size())
      tau[is].resize(3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector t = atom_list[is][ia]->position();
      tau[is][i++] = t.x;
      tau[is][i++] = t.y;
      tau[is][i++] = t.z;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_positions(const vector<vector<double> >& tau)
{
  assert(tau.size() == atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    assert(tau[is].size() == 3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      atom_list[is][ia]->set_position(
        D3vector(tau[is][i],tau[is][i+1],tau[is][i+2]));
      i += 3;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::get_velocities(vector<vector<double> >& vel) const
{
  if (vel.size() != atom_list.size())
    vel.resize(atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    if (vel[is].size() != 3*atom_list[is].size())
      vel[is].resize(3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector t = atom_list[is][ia]->velocity();
      vel[is][i++] = t.x;
      vel[is][i++] = t.y;
      vel[is][i++] = t.z;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_velocities(const vector<vector<double> >& vel)
{
  assert(vel.size() == atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    assert(vel[is].size() == 3*atom_list[is].size());
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      atom_list[is][ia]->set_velocity(
        D3vector(vel[is][i],vel[is][i+1],vel[is][i+2]));
      i += 3;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::reset_velocities(void)
{
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
      atom_list[is][ia]->set_velocity(D3vector(0.0, 0.0, 0.0));
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::rescale_velocities(double fac)
{
  vector<vector<double> > v;
  get_velocities(v);
  for ( int is = 0; is < v.size(); is++ )
  {
    for ( int ia = 0; ia < v[is].size(); ia++ )
      v[is][ia] *= fac;
  }
  set_velocities(v);
}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::vcm(void) const
{
  D3vector mvsum;
  double msum = 0.0;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double mass = species_list[is]->mass();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector v = atom_list[is][ia]->velocity();
      mvsum += mass * v;
      msum += mass;
    }
  }
  if ( msum == 0.0 ) return D3vector(0,0,0);
  return mvsum / msum;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::reset_vcm(void)
{
  D3vector vc = vcm();
  vector<vector<double> > v;
  get_velocities(v);
  // subtract center of mass velocity
  for ( int is = 0; is < v.size(); is++ )
  {
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      v[is][i++] -= vc.x;
      v[is][i++] -= vc.y;
      v[is][i++] -= vc.z;
    }
  }
  set_velocities(v);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::fold_in_ws(void)
{
  vector<vector<double> > p;
  get_positions(p);
  for ( int is = 0; is < p.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector pos(&p[is][3*ia]);
      cell_.fold_in_ws(pos);
      p[is][3*ia+0] = pos.x;
      p[is][3*ia+1] = pos.y;
      p[is][3*ia+2] = pos.z;
    }
  }
  set_positions(p);
}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::dipole(void) const
{
  D3vector sum;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double charge = species_list[is]->zval();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector p = atom_list[is][ia]->position();
      sum += charge * p;
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::sync()
{
#if USE_MPI
  // enforce consistency of positions and velocities on all tasks
  // broadcast positions and velocities of task 0 to all tasks
  vector<vector<double> > r,v;
  get_positions(r);
  get_velocities(v);
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    int m = r[is].size();
    double* p = &r[is][0];
    if ( ctxt_.onpe0() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    int m = v[is].size();
    double* p = &v[is][0];
    if ( ctxt_.onpe0() )
    {
      ctxt_.dbcast_send(m,1,p,m);
    }
    else
    {
      ctxt_.dbcast_recv(m,1,p,m,0,0);
    }
  }
  set_positions(r);
  set_velocities(v);

  // synchronize the unit cell
  if ( ctxt_.onpe0() )
  {
    double sbuf[9];
    for ( int i = 0; i < 9; i++ )
      sbuf[i] = cell_.amat(i);
    ctxt_.dbcast_send(9,1,sbuf,9);
  }
  else
  {
    double rbuf[9];
    ctxt_.dbcast_recv(9,1,rbuf,9,0,0);
    D3vector a0(rbuf[0],rbuf[1],rbuf[2]);
    D3vector a1(rbuf[3],rbuf[4],rbuf[5]);
    D3vector a2(rbuf[6],rbuf[7],rbuf[8]);
    // call UnitCell::set to recompute other UnitCell members
    cell_.set(a0,a1,a2);
  }
#endif
}
////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, const AtomSet &as )
{
  if ( as.context().onpe0() )
  {
    os << "<atomset>\n";
    os << as.cell();
    for ( int is = 0; is < as.species_list.size(); is++ )
    {
      os << *as.species_list[is];
    }
    for ( int is = 0; is < as.species_list.size(); is++ )
    {
      for ( int ia = 0; ia < as.atom_list[is].size(); ia++ )
      {
        os << *as.atom_list[is][ia];
      }
    }
  }
  os << "</atomset>\n";

  return os;
}
