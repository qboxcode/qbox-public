////////////////////////////////////////////////////////////////////////////////
//
// AtomSet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomSet.C,v 1.12 2006-03-07 07:02:28 fgygi Exp $

#include "AtomSet.h"
#include "Species.h"
#include "NameOf.h"
#include <iostream>
#include <algorithm>
#include <string>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
bool AtomSet::addSpecies(Species* sp, string name)
{
  if ( findSpecies(name) )
  {
    if ( ctxt_.onpe0() )
      cout << " AtomSet::addSpecies: species " << name
           << " is already defined" << endl;
    return false;
  }
  
  const double rcps = 1.0;
  sp->initialize(rcps);
 
  // create new entry in species list
  species_list.push_back(sp);
  spname.push_back(name);
  isp_[name] = species_list.size()-1;
  na_.insert(map<string,int>::value_type(name,0));
  atom_list.resize(atom_list.size()+1);
  is_[name] = spname.size()-1;
  
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
int AtomSet::size(void)
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
    int i = 0;
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      atom_list[is][ia]->set_velocity(D3vector(0.0, 0.0, 0.0));
      i += 3;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::sync()
{
#if USE_MPI
  // enforce consistency of positions on all tasks
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
#endif
}
////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, AtomSet &as )
{
  if ( as.context().onpe0() )
  {
    os << "<atomset>\n";
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
