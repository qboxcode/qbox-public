////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
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
// ExtForceSet.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "ExtForceSet.h"
#include "AtomicExtForce.h"
#include "PairExtForce.h"
#include "GlobalExtForce.h"
#include "Atom.h"
#include "AtomSet.h"
#include <iostream>
#include <iomanip>
#include <cstdlib> // atof
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ExtForceSet::~ExtForceSet(void)
{
  for ( int i = 0; i < extforce_list.size(); i++ )
    delete extforce_list[i];
}

////////////////////////////////////////////////////////////////////////////////
bool ExtForceSet::define_extforce(AtomSet &atoms, int argc, char **argv)
{
  enum extforce_type { unknown, atomic_type, pair_type, global_type }
    type = unknown;

  // argv[0] == "extforce"
  // argv[1] == "define"
  // argv[2] == {"atomic", "pair", "global"}
  // argv[3] == extforce_name
  // if extforce_type == atomic:
  //    argv[4] == {atom_name}
  //    argv[{5,6,7}] == {force vector}
  // else if extforce_type == pair:
  //    argv[{4,5}] == {atom names}
  //    argv[6] == {force_magnitude}
  // else if extforce_type == global:
  //    argv[{4,5,6}] == {force vector}

  if ( argc < 7 )
  {
    if ( MPIdata::onpe0() )
    {
      cout << " Use: extforce define atomic name atom fx fy fz"
           << endl;
      cout << " Use: extforce define pair name atom1 atom2 f"
           << endl;
      cout << " Use: extforce define global name fx fy fz"
           << endl;
    }
    return false;
  }
  const string extforce_type = argv[2];
  if ( extforce_type == "atomic" )
  {
    type = atomic_type;
  }
  else if ( extforce_type == "pair" )
  {
    type = pair_type;
  }
  else if ( extforce_type == "global" )
  {
    type = global_type;
  }
  else
  {
    if ( MPIdata::onpe0() )
      cout << " Incorrect extforce type " << extforce_type << endl;
    return false;
  }

  if ( type == atomic_type )
  {
    // extforce define atomic name atom fx fy fz

    if ( argc != 8 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for atomic extforce"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];

    double fx = atof(argv[5]);
    double fy = atof(argv[6]);
    double fz = atof(argv[7]);
    D3vector f(fx,fy,fz);

    Atom *a1 = atoms.findAtom(name1);

    if ( a1 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        cout << " ExtForceSet: could not find atom " << name1 << endl;
        cout << " ExtForceSet: could not define extforce" << endl;
      }
      return false;
    }

    // check if extforce is already defined
    bool found = false;
    ExtForce *pc = 0;
    for ( int i = 0; i < extforce_list.size(); i++ )
    {
      pc = extforce_list[i];
      assert(pc != 0);
      // check if an extforce with same name or with same atom is defined
      if ( pc->type() == "atomic" )
        found = ( pc->name() == name ) || ( pc->names(0) == name1 );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ExtForceSet: extforce is already defined:\n"
             << " cannot define extforce" << endl;
      return false;
    }
    else
    {
      AtomicExtForce *xf =
        new AtomicExtForce(name,name1,f);
      extforce_list.push_back(xf);
    }
  }
  else if ( type == pair_type )
  {
    // define pair name A B force_magnitude
    if ( argc != 7 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for pair extforce"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];

    // force magnitude
    double f = atof(argv[6]);

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);

    if ( a1 == 0 || a2 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        if ( a1 == 0 )
          cout << " ExtForceSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ExtForceSet: could not find atom " << name2 << endl;
        cout << " ExtForceSet: could not define extforce" << endl;
      }
      return false;
    }
    if ( name1 == name2 )
    {
      if ( MPIdata::onpe0() )
        cout << " ExtForceSet: cannot define pair extforce between "
             << name1 << " and " << name2 << endl;
      return false;
    }

    // check if extforce is already defined
    bool found = false;
    ExtForce *pc = 0;
    for ( int i = 0; i < extforce_list.size(); i++ )
    {
      pc = extforce_list[i];
      assert(pc != 0);
      // check if an extforce (name1,name2) or (name2,name1) is defined
      if ( pc->type() == "pair" )
        found = ( pc->names(0) == name1 && pc->names(1) == name2 ) ||
                ( pc->names(1) == name1 && pc->names(0) == name2 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ExtForceSet: extforce is already defined:\n"
             << " cannot define extforce" << endl;
      return false;
    }
    else
    {
      PairExtForce *xf =
        new PairExtForce(name,name1,name2,f);

      extforce_list.push_back(xf);
    }
  }
  else if ( type == global_type )
  {
    // define global name fx fy fz
    if ( argc != 7 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for global extforce"
             << endl;
      return false;
    }
    string name = argv[3];
    double fx = atof(argv[4]);
    double fy = atof(argv[5]);
    double fz = atof(argv[6]);
    D3vector f(fx,fy,fz);

    // check if extforce is already defined
    bool found = false;
    ExtForce *pc = 0;
    for ( int i = 0; i < extforce_list.size(); i++ )
    {
      pc = extforce_list[i];
      assert(pc != 0);
      // check if an extforce (name1,name2) or (name2,name1) is defined
      if ( pc->type() == "global" )
        found = ( pc->name() == name );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ExtForceSet: extforce is already defined:\n"
             << " cannot define extforce" << endl;
      return false;
    }
    else
    {
      GlobalExtForce *xf = new GlobalExtForce(name,f);
      extforce_list.push_back(xf);
    }
  }
  else
  {
    if ( MPIdata::onpe0() )
      cout << " ExtForceSet::set_constraint: internal error" << endl;
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ExtForceSet::set_extforce(int argc, char **argv)
{
  assert(argc==4||argc==6);
  // argv[0] == "extforce"
  // argv[1] == "set"
  // argv[2] == extforce_name
  // argv[3] == value
  // or argv[3,4,5] == value
  string name = argv[2];
  vector<double> value(argc-3);
  for ( int i = 0; i < argc-3; i++ )
    value[i] = atof(argv[i+3]);

  // check if constraint is defined and set its value
  bool found = false;
  vector<ExtForce*>::iterator i = extforce_list.begin();
  while ( !found && i != extforce_list.end() )
  {
    ExtForce *pc = *i;
    assert(pc != 0);

    if ( pc->name() == name )
    {
      found = true;
      pc->set_value(value);
    }
    i++;
  }

  if ( !found )
  {
    if ( MPIdata::onpe0() )
      cout << " ExtForceSet: no such extforce" << endl;
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ExtForceSet::delete_extforce(int argc, char **argv)
{
  assert(argc==3);
  // argv[0] == "extforce"
  // argv[1] == "delete"
  // argv[2] == "extforce_name"
  string name = argv[2];

  bool found = false;
  // note next loop in reverse: avoid use of invalidated iterators
  // after erase operation

  vector<ExtForce*>::iterator i = extforce_list.begin();
  while ( !found && i != extforce_list.end() )
  {
    ExtForce *pc = *i;
    assert(pc != 0);

    // note structure of if else test to avoid incrementing
    // invalidated iterator after erase (see Meyers STL, p.45)
    if ( pc->name() == name )
    {
      found = true;
      delete pc;

      // remove pointer from the list
      // note: iterator is incremented before erasing, remains valid
      extforce_list.erase(i++);
    }
    else
    {
      i++;
    }
  }

  if ( !found )
  {
    if ( MPIdata::onpe0() ) cout << " No such extforce" << endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void ExtForceSet::list_extforces(ostream &os)
{
  if ( !extforce_list.empty() )
  {
    const double au2pN = 82388.5845379; // Ha/Bohr -> pN conversion
    D3vector sum;
    os << "<extforce_set>" << endl;
    for ( int i = 0; i < extforce_list.size(); i++ )
    {
      ExtForce *x = extforce_list[i];
      os << *x << endl;
      sum += x->sum_contrib();
      const double mag = x->magnitude();
      os << "  <!-- magnitude: " << mag << " (a.u.) = "
         << au2pN*mag << " pN -->" << endl;
    }
    os << "  <!-- sum: " << sum << " -->" << endl;
    os << "  <!-- sum magnitude: " << length(sum) << " (a.u.) = "
       << au2pN*length(sum) << " pN -->" << endl;
    os << "</extforce_set>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
double ExtForceSet::energy(const vector<vector<double> > &r,
                           vector<vector<double> > &f) const
{
  for ( int is = 0; is < f.size(); is++ )
    for ( int i = 0; i < f[is].size(); i++ )
      f[is][i] = 0.0;
  double sum = 0;
  for ( int i = 0; i < extforce_list.size(); i++ )
  {
    sum += extforce_list[i]->energy(r,f);
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void ExtForceSet::setup(AtomSet& atoms)
{
  for ( int i = 0; i < extforce_list.size(); i++ )
  {
    extforce_list[i]->setup(atoms);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExtForceSet::reset(void)
{
  for ( int i = 0; i < extforce_list.size(); i++ )
    delete extforce_list[i];
  extforce_list.resize(0);
}
