////////////////////////////////////////////////////////////////////////////////
//
// ConstraintSet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintSet.C,v 1.1 2005-06-27 22:34:46 fgygi Exp $

#include "ConstraintSet.h"
#include "DistanceConstraint.h"
#include "AngleConstraint.h"
#include "TorsionConstraint.h"
//#include "MultiDistanceConstraint.h"
#include "Atom.h"
#include "AtomSet.h"
#include "Context.h"
#include <iostream>
#include <iomanip>
using namespace std;

const int constraints_maxiter = 10;

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::set_constraint(AtomSet &atoms, int argc, char **argv)
{
  enum constraint_type { unknown, distance_type, multidistance_type, 
                         angle_type, torsion_type }
    type = unknown;
  const double distance_tolerance = 1.0e-7;
  const double angle_tolerance = 1.0e-4;
  const bool onpe0 = ctxt_.onpe0();

  if ( argc < 2 )
  {
    if ( onpe0 )
    {
      cout << " use: constraint distance name1 name2 distance [velocity]" 
           << endl;
      cout << "      constraint angle name1 name2 name3 angle [velocity]"
           << endl;
      cout << "      constraint torsion name1 name2 name3 name4 angle"
           << " [velocity]" 
           << endl;
      //cout << "      constraint multidistance alpha_i name1_i name2_i"
      //     << " distance [velocity]" 
      //<< endl;
    }
    return 1;
  }
    
  if ( !strcmp(argv[1],"distance") )
  {
    type = distance_type;
  }
  else if ( !strcmp(argv[1],"multidistance") )
  { 
    type = multidistance_type;
    if ( onpe0 ) cout << " multidistance  constraint not implemented" << endl;
    return false;
  }
  else if ( !strcmp(argv[1],"angle") )
  {
    type = angle_type;
  }
  else if ( !strcmp(argv[1],"torsion") )
  {
    type = torsion_type;
  }
  else
  {
    if ( onpe0 ) cout << " Incorrect constraint type" << endl;
    return false;
  }
 
  if ( type == distance_type )
  {
    // constraint distance A B dist
    // constraint distance A B dist velocity

    if ( argc < 5 || argc > 6 )
    {
      if ( onpe0 ) 
        cout << " Incorrect number of arguments for distance constraint" 
             << endl;
      return false;
    }
    double distance, velocity=0.0, distance_tolerance=1.e-7;
    string name1 = argv[2];
    string name2 = argv[3];
 
    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
 
    if ( a1 == 0 || a2 == 0 )
    {
      if ( a1 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name1 << endl;
      if ( a2 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name2 << endl;
 
      if ( onpe0 ) 
        cout << " ConstraintSet: could not define constraint " << endl;
      return 1;
    }
    if ( name1 == name2 )
    {
      if ( onpe0 ) 
        cout << " ConstraintSet: cannot define distance constraint between "
             << name1 << " and " << name2 << endl;
      return 1;
    }
 
    distance = atof(argv[4]);
    if ( argc == 6 )
    {
      velocity = atof(argv[5]);
    }
 
    if ( distance <= 0.0 )
    {
      if ( onpe0 )
        cout << " ConstraintSet: distance must be positive " << endl
             << " ConstraintSet: could not define constraint " << endl;
      return 1;
    }
  
    // check if equivalent constraint is already defined
    bool found = false;
    Constraint *pc = 0;
    for ( int i = 0; i < constraint_list.size(); i++ )
    {
      pc = constraint_list[i];
      assert(pc != 0);
      // check if a constraint (name1,name2) or (name2,name1) is defined
      if ( pc->type() == "distance" &&
           ( pc->names(0) == name1 && pc->names(1) == name2 ) ||
           ( pc->names(1) == name1 && pc->names(0) == name2 ) )
         found = true;
    }
 
    if ( found )
    {
      if ( onpe0 )
        cout << " <!-- ConstraintSet:set_constraint: distance constraint "
             << name1 << " " << name2
             << " was found" << endl
             << " ConstraintSet: modifying distance constraint -->" << endl;
      pc->set_value(distance);
      pc->set_velocity(velocity);
    }
    else
    {
      DistanceConstraint *c =
        new DistanceConstraint(name1,name2,distance,
                               velocity,distance_tolerance);
    
      constraint_list.push_back(c);
    }
  }
  else if ( type == angle_type )
  {
    // constraint angle A B C angle
    // constraint angle A B C angle velocity

    if ( argc < 6  || argc > 7 )
    {
      if ( onpe0 )
        cout << " Incorrect number of arguments for angle constraint"
             << endl;
      return false;
    }
    string name1 = argv[2];
    string name2 = argv[3];
    string name3 = argv[4];
 
    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);
 
    if ( a1 == 0 || a2 == 0 || a3 == 0 )
    {
      if ( a1 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name1 << endl;
      if ( a2 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name2 << endl;
      if ( a3 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name3 << endl;
 
      if ( onpe0 ) 
        cout << " ConstraintSet: could not define constraint " << endl;
      return 1;
    }
    if ( name1 == name2 || name1 == name3 || name2 == name3)
    {
      if ( onpe0 ) 
        cout << " ConstraintSet: cannot define angle constraint between "
             << name1 << " " << name2 << " and " << name3 << endl;
      return 1;
    }
 
    const double angle = atof(argv[5]);
    double velocity;
    if ( argc == 7 )
    {
      velocity = atof(argv[6]);
    }
 
    if ( angle < 0.0 || angle > 180.0 )
    {
      if ( onpe0 )
        cout << " ConstraintSet: angle must be in [0,180]" << endl
             << " ConstraintSet: could not define constraint " << endl;
      return 1;
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
      if ( pc->type() == "angle" &&
           ( pc->names(0) == name1 && 
             pc->names(1) == name2 && 
             pc->names(2) == name3) ||
           ( pc->names(0) == name3 && 
             pc->names(1) == name2 &&
             pc->names(2) == name1) )
         found = true;
    }
 
    if ( found )
    {
      if ( onpe0 )
        cout << " <!-- ConstraintSet:set_constraint: an angle constraint "
             << name1 << " " << name2 << " " << name3
             << " was found" << endl
             << " ConstraintSet: modifying angle constraint -->" << endl;
        pc->set_value(angle);
        pc->set_velocity(velocity);
    }
    else
    {
      AngleConstraint *c =
      new AngleConstraint(name1,name2,name3,angle,velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
  else if ( type == torsion_type )
  {
    // constraint torsion A B C D angle
    // constraint torsion A B C D angle velocity

    if ( argc < 7  || argc > 8 )
    {
      if ( onpe0 )
        cout << " Incorrect number of arguments for angle constraint"
             << endl;
      return false;
    }
    string name1 = argv[2];
    string name2 = argv[3];
    string name3 = argv[4];
    string name4 = argv[5];
 
    Atom *a1 = atoms.findAtom(name1);
    Atom *a2 = atoms.findAtom(name2);
    Atom *a3 = atoms.findAtom(name3);
    Atom *a4 = atoms.findAtom(name4);
 
    if ( a1 == 0 || a2 == 0 || a3 == 0 || a4 == 0 )
    {
      if ( a1 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name1 << endl;
      if ( a2 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name2 << endl;
      if ( a3 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name3 << endl;
      if ( a4 == 0 )
        if ( onpe0 ) 
          cout << " ConstraintSet: could not find atom " << name4 << endl;
 
      if ( onpe0 ) 
        cout << " ConstraintSet: could not define constraint " << endl;
      return 1;
    }
    if ( name1 == name2 || name1 == name3 || name1 == name4 ||
         name2 == name3 || name2 == name4 || name3 == name4 )
    {
      if ( onpe0 ) 
        cout << " ConstraintSet: cannot define torsion constraint using "
             << name1 << " " << name2 << " " << name3 << " " << name4 << endl;
      return 1;
    }
 
    double angle = atof(argv[6]);
    if ( angle > 180.0 )
      while ( angle > 180.0 ) angle -= 360.0;
    else if ( angle < -180.0 )
      while ( angle < -180.0 ) angle += 360.0;

    double velocity;
    if ( argc == 8 )
    {
      velocity = atof(argv[7]);
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
      if ( pc->type() == "angle" &&
           ( pc->names(0) == name1 && 
             pc->names(1) == name2 && 
             pc->names(2) == name3 && 
             pc->names(3) == name4) ||
           ( pc->names(0) == name4 && 
             pc->names(1) == name3 &&
             pc->names(2) == name2 &&
             pc->names(3) == name1) )
         found = true;
    }
 
    if ( found )
    {
      if ( onpe0 )
        cout << " <!-- ConstraintSet:add_constraint: a torsion constraint "
             << name1 << " " << name2 << " " << name3 << " " << name4
             << " is already defined" << endl
             << " ConstraintSet: modifying torsion constraint -->" << endl;
      pc->set_value(angle);
      pc->set_velocity(velocity);
    }
    else
    {
      TorsionConstraint *c =
      new TorsionConstraint(name1,name2,name3,name4,
                            angle,velocity,angle_tolerance);
      constraint_list.push_back(c);
    }
  }
#if 0
  else if ( type == multidistance_type )
  {
    // constraint multidistance alpha1 A1 B1 alpha2 A2 B2 ... d 
    // constraint multidistance alpha1 A1 B1 alpha2 A2 B2 ... d v

    vector<string> m_name1, m_name2;
    vector<double> m_alpha;
    int nc;
    if ( argc < 6 )
    {
      if ( onpe0 )
        cout << " Incorrect number of arguments for multidistance constraint"
             << endl;
      return false;
    }
    if ( argc % 3 == 0 )   // no velocity term included
    {
      distance = atof(argv[argc-1]);
      nc = argc-1;
    }
    else if ( (argc-1) % 3 == 0 )  // velocity term included
    {
      distance = atof(argv[argc-2]);
      velocity = atof(argv[argc-1]);
      nc = argc-2;
    }
    else
    {
      if ( onpe0 )
        cout << " Incorrect number of arguments for multidistance constraint"
             << endl;
      return false;
    }

    for ( int ic=2; ic < nc; ic=ic+3 )
    {
      double alpha12 = atof(argv[ic]);   
      name1 = argv[ic+1];
      name2 = argv[ic+2];

      Atom *a1 = atoms.findAtom(name1);
      Atom *a2 = atoms.findAtom(name2);

      if ( a1 == 0 || a2 == 0 )
      {
        if ( a1 == 0 )
          if ( onpe0 )
            cout << " ConstraintSet: could not find atom " << name1 << endl;
        if ( a2 == 0 )
          if ( onpe0 )
            cout << " ConstraintSet: could not find atom " << name2 << endl;

        if ( onpe0 )
          cout << " ConstraintSet: could not define constraint " << endl;
        return 1;
      }
      if ( name1 == name2 )
      {
        if ( onpe0 )
          cout << " ConstraintSet: cannot define distance constraint between "
               << name1 << " and " << name2 << endl;
        return 1;
      }

     // put the atom names in alphabetical order

      if ( name1 < name2 ) 
      {
        m_name1.push_back(name1);
        m_name2.push_back(name2);
      }
      else
      {
        m_name1.push_back(name2);
        m_name2.push_back(name1);
      }
      m_alpha.push_back(alpha12);
    }

    MultiDistanceConstraint *c =
      new MultiDistanceConstraint(m_alpha,m_name1,m_name2,distance,
                                  velocity,tolerance);

    constraint_list.push_back(c);

  }
#endif
  else
  {
    if ( onpe0 )
      cout << "ConstraintSet::add_constraint: internal error" << endl;
    return false;
  }
 
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool ConstraintSet::del_constraint(int argc, char **argv)
{
  enum constraint_type { unknown, distance_type, multidistance_type, 
                         angle_type, torsion_type }
    type = unknown;
  string name1, name2, name3, name4;
  const bool onpe0 = ctxt_.onpe0();

  if ( argc < 2 )
  {
    if ( onpe0 )
    {
      cout << " use: del_constraint constraint_type atom_name [atom_name]" 
           << endl;
      cout << "      del_constraint all" << endl;
    }
    return 1;
  }
    
  if ( !strcmp(argv[1],"all") )
  {
    for ( int i = 0; i < constraint_list.size(); i++)
    {
      delete constraint_list[i];  
    }
    constraint_list.resize(0);
    return true;
  }
  else if ( !strcmp(argv[1],"distance") )
  {
    type = distance_type;
  }
  else if ( !strcmp(argv[1],"multidistance") )
  {
    type = multidistance_type;
    if ( onpe0 ) cout << 
         " multidistance constraint deletion not implemented" << endl;
    return false;
  }
  else if ( !strcmp(argv[1],"angle") )
  {
    type = angle_type;
  }
  else if ( !strcmp(argv[1],"torsion") )
  {
    type = torsion_type;
  }
  else
  {
    if ( onpe0 ) cout << " Incorrect constraint type" << endl;
    return false;
  }
 
  if ( type == distance_type )
  {
    // del_constraint distance X Y
    if ( argc != 4 )
    {
      if ( onpe0 ) 
        cout << " Incorrect number of arguments" << endl;
      return false;
    }
    name1 = argv[2];
    name2 = argv[3];
 
    bool found = false;
    // note next loop in reverse: avoid use of invalidated iterators
    // after erase operation
    
    vector<Constraint*>::iterator i = constraint_list.begin();
    while ( !found && i != constraint_list.end() )
    {
      Constraint *pc = *i;
      assert(pc != 0);
      // check if a constraint (name1,name2) or (name2,name1) is defined
      
      // note structure of if else test to avoid incrementing
      // invalidated iterator after erase (see Meyers STL, p.45)
      if ( pc->type() == "distance" &&
           ( pc->names(0) == name1 && pc->names(1) == name2 ) ||
           ( pc->names(1) == name1 && pc->names(0) == name2 ) )
      {
        found = true;
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
      if ( onpe0 ) cout << " No such constraint" << endl;
  }
    
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::list_constraints(ostream &os)
{
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    Constraint *c = constraint_list[i];
    os << *c << endl;
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
  const bool onpe0 = ctxt_.onpe0();
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
  
  for ( int i = 0; i < constraint_list.size(); i++ )  
  {                                                   
    Constraint *c = constraint_list[i];               
    if ( c->velocity() != 0.0 )                       
    {                                                 
      if ( onpe0 )                                    
        cout << *c << endl;                           
    }                                                 
  }
  
  if ( !done )
  {
    if ( onpe0 )
      cout << " <!-- ConstraintSet: could not enforce position constraints in "
           << constraints_maxiter << " iterations -->" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ConstraintSet::enforce_v(const vector<vector<double> > &r0,
                              vector<vector<double> > &v0)
{
  const bool onpe0 = ctxt_.onpe0();
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
    if ( onpe0 )
      cout << " <!-- ConstraintSet: could not enforce velocity constraints in "
           << constraints_maxiter << " iterations -->" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
double ConstraintSet::projection(const vector<vector<double> > &r0,
 const vector<vector<double> > &x)
{
  double sum = 0.0;
  for ( int i = 0; i < constraint_list.size(); i++ )
  {
    sum += constraint_list[i]->projection(r0,x);
  }
  return sum;
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
