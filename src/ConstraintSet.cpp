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
// ConstraintSet.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "ConstraintSet.h"
#include "PositionConstraint.h"
#include "DistanceConstraint.h"
#include "AngleConstraint.h"
#include "TorsionConstraint.h"
#include "Atom.h"
#include "AtomSet.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cstdlib> // atof
#include <cstring> // strcmp
using namespace std;

const int constraints_maxiter = 50;

////////////////////////////////////////////////////////////////////////////////
ConstraintSet::~ConstraintSet(void)
{
  for ( int ic = 0; ic < constraint_list.size(); ic++ )
    delete constraint_list[ic];
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::define_constraint(AtomSet &atoms, int argc, char **argv)
{
  enum constraint_type { unknown, position_type, distance_type,
                         angle_type, torsion_type }
    type = unknown;
  const double position_tolerance = 1.0e-7;
  const double distance_tolerance = 1.0e-7;
  const double angle_tolerance = 1.0e-4;

  // argv[0] == "constraint"
  // argv[1] == "define"
  // argv[2] == {"position", "distance", "angle", "torsion"}
  // argv[3] == constraint name
  // argv[4-(5,6,7)] == atom names
  // argv[{5,6,7}] == {distance,angle,angle}
  // argv[{6,7,8}] == velocity

  if ( argc < 2 )
  {
    if ( MPIdata::onpe0() )
    {
      cout << " Use: constraint define position constraint_name atom_name"
           << endl;
      cout << " Use: constraint define distance constraint_name "
           << "atom_name1 atom_name2 distance_value [velocity]"
           << endl;
      cout << "      constraint define angle constraint_name "
           << "name1 name2 name3 angle_value [velocity]"
           << endl;
      cout << "      constraint define torsion constraint_name "
           << "name1 name2 name3 name4 angle_value"
           << " [velocity] "
           << endl;
    }
    return false;
  }
  const string constraint_type = argv[2];
  if ( constraint_type == "position" )
  {
    type = position_type;
  }
  else if ( constraint_type == "distance" )
  {
    type = distance_type;
  }
  else if ( constraint_type == "angle" )
  {
    type = angle_type;
  }
  else if ( constraint_type == "torsion" )
  {
    type = torsion_type;
  }
  else
  {
    if ( MPIdata::onpe0() )
      cout << " Incorrect constraint type " << constraint_type << endl;
    return false;
  }

  if ( type == position_type )
  {
    // define position name A

    if ( argc != 5 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for position constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];

    Atom *a1 = atoms.findAtom(name1);

    if ( a1 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        cout << " ConstraintSet: could not find atom " << name1 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }

    // check if constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint with same name or with same atom is defined
      if ( pc->type() == "position" )
        found = ( pc->name() == name ) || ( pc->names(0) == name1 );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: constraint is already defined:\n"
             << " cannot define constraint" << endl;
      return false;
    }
    else
    {
      PositionConstraint *c =
        new PositionConstraint(name,name1,position_tolerance);
      constraint_list.push_back(c);
    }
  }
  else if ( type == distance_type )
  {
    // define distance name A B value
    // define distance name A B value velocity

    if ( argc < 7 || argc > 8 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for distance constraint"
             << endl;
      return false;
    }
    double distance, velocity=0.0;
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);

    if ( a1 == 0 || a2 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }
    if ( name1 == name2 )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: cannot define distance constraint between "
             << name1 << " and " << name2 << endl;
      return false;
    }

    // define distance value
    if ( !strcmp(argv[6],"*") )
    {
      // use current distance
      distance = length(a1->position()-a2->position());
      if ( MPIdata::onpe0() )
        cout << "ConstraintSet::define_constraint: using current distance "
             << distance << endl;
    }
    else
      distance = atof(argv[6]);

    if ( argc == 8 )
    {
      velocity = atof(argv[7]);
    }

    if ( distance <= 0.0 )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: distance must be positive" << endl
             << " ConstraintSet: could not define constraint" << endl;
      return false;
    }

    // check if constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint (name1,name2) or (name2,name1) is defined
      if ( pc->type() == "distance" )
        found = ( pc->names(0) == name1 && pc->names(1) == name2 ) ||
                ( pc->names(1) == name1 && pc->names(0) == name2 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: a distance constraint named " << name << endl
             << " or involving atoms " << name1 << " " << name2 << endl
             << " is already defined" << endl
             << " ConstraintSet: cannot define constraint" << endl;
      return false;
    }
    else
    {
      DistanceConstraint *c =
        new DistanceConstraint(name,name1,name2,distance,
                               velocity,distance_tolerance);

      constraint_list.push_back(c);
    }
  }
  else if ( type == angle_type )
  {
    // constraint define angle name A B C value
    // constraint define angle name A B C value velocity

    if ( argc < 8  || argc > 9 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for angle constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];
    string name3 = argv[6];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);

    if ( a1 == 0 || a2 == 0 || a3 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        if ( a3 == 0 )
          cout << " ConstraintSet: could not find atom " << name3 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }

    if ( name1 == name2 || name1 == name3 || name2 == name3)
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: cannot define angle constraint between "
             << name1 << " " << name2 << " and " << name3 << endl;
      return false;
    }

    // define angle value
    double angle = 0.0;
    if ( !strcmp(argv[7],"*") )
    {
      // use current angle
      D3vector r12(a1->position()-a2->position());
      D3vector r32(a3->position()-a2->position());
      if ( norm2(r12) == 0.0 || norm2(r32) == 0.0 )
      {
        if ( MPIdata::onpe0() )
        {
          cout << " ConstraintSet: cannot define angle constraint." << endl;
          cout << " ConstraintSet: atoms are too close" << endl;
          return false;
        }
      }
      const double sp = normalized(r12) * normalized(r32);
      const double c = max(-1.0,min(1.0,sp));
      angle = (180.0/M_PI)*acos(c);

      if ( MPIdata::onpe0() )
        cout << "ConstraintSet::define_constraint: using current angle "
             << angle << endl;
    }
    else
      angle = atof(argv[7]);

    double velocity = 0.0;
    if ( argc == 9 )
    {
      velocity = atof(argv[8]);
    }

    if ( angle < 0.0 || angle > 180.0 )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: angle must be in [0,180]" << endl
             << " ConstraintSet: could not define constraint" << endl;
      return false;
    }

    // check if equivalent constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint (name1,name2,name3) or
      // (name3,name2,name1) is defined
      if ( pc->type() == "angle" )
        found = ( pc->names(0) == name1 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name3 ) ||
                ( pc->names(0) == name3 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name1 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: an angle constraint named " << name << endl
             << " or involving atoms "
             << name1 << " " << name2 << " " << name3 << endl
             << " is already defined" << endl
             << " ConstraintSet: cannot define constraint" << endl;
      return false;
    }
    else
    {
      AngleConstraint *c =
      new AngleConstraint(name, name1,name2,name3,angle,
        velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
  else if ( type == torsion_type )
  {
    // constraint define torsion name A B C D angle
    // constraint define torsion name A B C D angle velocity

    if ( argc < 9  || argc > 10 )
    {
      if ( MPIdata::onpe0() )
        cout << " Incorrect number of arguments for torsion constraint"
             << endl;
      return false;
    }
    string name = argv[3];
    string name1 = argv[4];
    string name2 = argv[5];
    string name3 = argv[6];
    string name4 = argv[7];

    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);
    Atom *a4 = atoms.findAtom(name4);

    if ( a1 == 0 || a2 == 0 || a3 == 0 || a4 == 0 )
    {
      if ( MPIdata::onpe0() )
      {
        if ( a1 == 0 )
          cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          cout << " ConstraintSet: could not find atom " << name2 << endl;
        if ( a3 == 0 )
          cout << " ConstraintSet: could not find atom " << name3 << endl;
        if ( a4 == 0 )
          cout << " ConstraintSet: could not find atom " << name4 << endl;
        cout << " ConstraintSet: could not define constraint" << endl;
      }
      return false;
    }
    if ( name1 == name2 || name1 == name3 || name1 == name4 ||
         name2 == name3 || name2 == name4 || name3 == name4 )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: cannot define torsion constraint using "
             << name1 << " " << name2 << " " << name3 << " " << name4
             << endl;
      return false;
    }

    // define angle value
    double angle = 0.0;
    if ( !strcmp(argv[8],"*") )
    {
      // use current angle
      D3vector r12(a1->position()-a2->position());
      D3vector r32(a3->position()-a2->position());
      D3vector r43(a4->position()-a3->position());
      if ( norm2(r12) == 0.0 || norm2(r32) == 0.0 || norm2(r43) == 0.0 )
      {
        if ( MPIdata::onpe0() )
        {
          cout << " ConstraintSet: cannot define torsion constraint." << endl;
          cout << " ConstraintSet: atoms are too close" << endl;
          return false;
        }
      }

      D3vector e12(normalized(r12));
      D3vector e32(normalized(r32));
      D3vector e23(-e32);
      D3vector e43(normalized(r43));

      const double sin123 = length(e12^e32);
      const double sin234 = length(e23^e43);

      if ( sin123 != 0.0 && sin234 != 0.0 )
      {
        D3vector e123 = normalized(e12^e32);
        D3vector e234 = normalized(e23^e43);
        double cc = max(min(e123*e234,1.0),-1.0);
        double ss = max(min((e123^e234)*e32,1.0),-1.0);
        angle = (180.0/M_PI) * atan2(ss,cc);
      }

      if ( MPIdata::onpe0() )
        cout << "ConstraintSet::define_constraint: using current angle "
             << angle << endl;
    }
    else
      angle = atof(argv[8]);

    if ( angle > 180.0 )
      while ( angle > 180.0 ) angle -= 360.0;
    else if ( angle < -180.0 )
      while ( angle < -180.0 ) angle += 360.0;

    double velocity = 0.0;
    if ( argc == 10 )
    {
      velocity = atof(argv[9]);
    }

    // check if equivalent constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if an equivalent constraint (name1,name2,name3,name4) or
      // (name4,name3,name2,name1) is defined
      if ( pc->type() == "torsion" )
        found = ( pc->names(0) == name1 &&
                  pc->names(1) == name2 &&
                  pc->names(2) == name3 &&
                  pc->names(3) == name4 ) ||
                ( pc->names(0) == name4 &&
                  pc->names(1) == name3 &&
                  pc->names(2) == name2 &&
                  pc->names(3) == name1 ) ||
                ( pc->name() == name );
    }

    if ( found )
    {
      if ( MPIdata::onpe0() )
        cout << " ConstraintSet: a torsion constraint named " << name << endl
             << " or involving atoms "
             << name1 << " " << name2 << " " << name3 << " " << name4 << endl
             << " is already defined" << endl
             << " ConstraintSet: cannot define constraint" << endl;
      return false;
    }
    else
    {
      TorsionConstraint *c =
      new TorsionConstraint(name,name1,name2,name3,name4,
                            angle,velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
  else
  {
    if ( MPIdata::onpe0() )
      cout << " ConstraintSet::set_constraint: internal error" << endl;
    return false;
  }

  // update total number of blocked degrees of freedom
  ndofs_ += constraint_list.back()->dofs();

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::set_constraint(int argc, char **argv)
{
  assert(argc==4||argc==5);
  // argv[0] == "constraint"
  // argv[1] == "set"
  // argv[2] == constraint_name
  // argv[3] == value
  // argv[4] (optional) == velocity
  string name = argv[2];
  const double value = atof(argv[3]);
  double velocity = 0.0;
  const bool set_velocity = ( argc == 5 );
  if ( set_velocity ) velocity = atof(argv[4]);

    // check if constraint is already defined
  bool found = false;
  vector<Constraint*>::iterator i = constraint_list.begin();
  while ( !found && i != constraint_list.end() )
  {
    Constraint *pc = *i;
    assert(pc != 0);

    if ( pc->name() == name )
    {
      found = true;
      pc->set_value(value);
      if ( set_velocity ) pc->set_velocity(velocity);
    }
    i++;
  }

  if ( !found )
  {
    if ( MPIdata::onpe0() )
      cout << " ConstraintSet: no such constraint" << endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::delete_constraint(int argc, char **argv)
{
  assert(argc==3);
  // argv[0] == "constraint"
  // argv[1] == "delete"
  // argv[2] == constraint_name
  string name = argv[2];

  bool found = false;
  // note next loop in reverse: avoid use of invalidated iterators
  // after erase operation

  vector<Constraint*>::iterator i = constraint_list.begin();
  while ( !found && i != constraint_list.end() )
  {
    Constraint *pc = *i;
    assert(pc != 0);

    // note structure of if else test to avoid incrementing
    // invalidated iterator after erase (see Meyers STL, p.45)
    if ( pc->name() == name )
    {
      found = true;

      // update total number of blocked degrees of freedom
      ndofs_ -= pc->dofs();

      delete pc;

      // remove constraint pointer from the list
      // note: iterator is incremented before erasing, remains valid
      constraint_list.erase(i++);
    }
    else
    {
      i++;
    }
  }

  if ( !found )
  {
    if ( MPIdata::onpe0() ) cout << " No such constraint" << endl;
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::list_constraints(ostream &os)
{
  if ( !constraint_list.empty() )
  {
    os << " <constraint_set>" << endl;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      Constraint *c = constraint_list[i];
      os << *c << endl;
    }
    os << " </constraint_set>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce(AtomSet& atoms)
{
  vector<vector<double> > r0,rp,v0;
  setup(atoms);
  atoms.get_positions(r0);
  rp=r0;
  atoms.get_velocities(v0);
  enforce_r(r0,rp);
  atoms.set_positions(rp);
  enforce_v(r0,v0);
  atoms.set_velocities(v0);
}
////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce_r(const vector<vector<double> > &r0,
                              vector<vector<double> > &rp)
{
  int iter = 0;
  bool done = false;
  while ( !done && (iter < constraints_maxiter) )
  {
    done = true;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      Constraint *c = constraint_list[i];
      bool b = c->enforce_r(r0,rp);
      done &= b;
    }
    iter++;
  }

  if ( !done )
  {
    if ( MPIdata::onpe0() )
      cout << " ConstraintSet: could not enforce position constraints in "
           << constraints_maxiter << " iterations" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce_v(const vector<vector<double> > &r0,
                              vector<vector<double> > &v0)
{
  int iter = 0;
  bool done = false;
  while ( !done && (iter < constraints_maxiter) )
  {
    done = true;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      bool b = constraint_list[i]->enforce_v(r0,v0);
      done &= b;
    }
    iter++;
  }

  if ( !done )
  {
    if ( MPIdata::onpe0() )
      cout << " ConstraintSet: could not enforce velocity constraints in "
           << constraints_maxiter << " iterations" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::compute_forces(const vector<vector<double> > &r0,
 const vector<vector<double> > &f)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->compute_force(r0,f);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::update_constraints(double dt)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->update(dt);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::setup(AtomSet& atoms)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    constraint_list[i]->setup(atoms);
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::reset(void)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
    delete constraint_list[i];
  ndofs_ = 0;
  constraint_list.resize(0);
}
