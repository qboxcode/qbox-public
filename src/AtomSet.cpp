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
// AtomSet.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "AtomSet.h"
#include "Species.h"
#include "NameOf.h"
#include "MPIdata.h"
#include "sampling.h"
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
    if ( MPIdata::onpe0() )
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

  if ( MPIdata::onpe0() )
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
    if ( MPIdata::onpe0() )
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
    if ( MPIdata::onpe0() )
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
  if ( MPIdata::onpe0() )
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
      if ( MPIdata::onpe0() )
        cout << *atom_list[is][ia];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::listSpecies(void) const
{
  for ( int is = 0; is < species_list.size(); is++ )
  {
    if ( MPIdata::onpe0() )
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
void AtomSet::sync_positions(vector<vector<double> >& tau)
{
  for ( int is = 0; is < tau.size(); is++ )
  {
    int m = tau[is].size();
    double* p = &tau[is][0];
    MPI_Bcast(p,m,MPI_DOUBLE,0,MPIdata::comm());
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_positions(vector<vector<double> >& tau)
{
  sync_positions(tau);
  assert(tau.size() == atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    assert(tau[is].size() == 3*atom_list[is].size());
    for ( int ia = 0, i = 0; ia < atom_list[is].size(); ia++ )
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
void AtomSet::sync_velocities(vector<vector<double> >& vel)
{
  for ( int is = 0; is < vel.size(); is++ )
  {
    int m = vel[is].size();
    double* p = &vel[is][0];
    MPI_Bcast(p,m,MPI_DOUBLE,0,MPIdata::comm());
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_velocities(vector<vector<double> >& vel)
{
  sync_velocities(vel);
  assert(vel.size() == atom_list.size());
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    assert(vel[is].size() == 3*atom_list[is].size());
    for ( int ia = 0, i = 0; ia < atom_list[is].size(); ia++ )
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
    for ( int i = 0; i < v[is].size(); i++ )
      v[is][i] *= fac;
  }
  set_velocities(v);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::randomize_positions(double amplitude)
{
  // add random displacements to positions using
  // random numbers from a normal distribution scaled
  // by the amplitude parameter
  vector<vector<double> > r;
  get_positions(r);
  for ( int is = 0; is < r.size(); is++ )
  {
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      // draw pairs of unit variance gaussian random variables
      double xi0, xi1, xi2, xi3; // xi3 not used
      normal_dev(&xi0,&xi1);
      normal_dev(&xi2,&xi3);
      r[is][3*ia+0] += amplitude * xi0;
      r[is][3*ia+1] += amplitude * xi1;
      r[is][3*ia+2] += amplitude * xi2;
    }
  }
  set_positions(r);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::randomize_velocities(double temp)
{
  // initialize velocities with random numbers from a Maxwell-Boltzmann
  // distribution at temperature temp
  const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
  vector<vector<double> > v;
  get_velocities(v);
  for ( int is = 0; is < v.size(); is++ )
  {
    const double m = species_list[is]->mass() * 1822.89;
    const double width = sqrt( boltz * temp / m );
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      // draw pairs of unit variance gaussian random variables
      double xi0, xi1, xi2, xi3; // xi3 not used
      normal_dev(&xi0,&xi1);
      normal_dev(&xi2,&xi3);
      v[is][3*ia+0] = width * xi0;
      v[is][3*ia+1] = width * xi1;
      v[is][3*ia+2] = width * xi2;
    }
  }
  set_velocities(v);
  reset_vcm();
}

////////////////////////////////////////////////////////////////////////////////
D3vector AtomSet::rcm(void) const
{
  D3vector mrsum;
  double msum = 0.0;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double mass = species_list[is]->mass();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector r = atom_list[is][ia]->position();
      mrsum += mass * r;
      msum += mass;
    }
  }
  if ( msum == 0.0 ) return D3vector(0,0,0);
  return mrsum / msum;
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
void AtomSet::reset_rotation(void)
{
  // check for special case of zero or one atom
  if ( size() < 2 ) return;

  D3vector rc = rcm();
  D3vector vc = vcm();
  vector<vector<double> > rt;
  get_positions(rt);
  vector<vector<double> > vt;
  get_velocities(vt);

  // compute angular momentum w.r.t. the center of mass
  D3vector L;
  for ( int is = 0; is < vt.size(); is++ )
  {
    double mass = species_list[is]->mass();
    for ( int ia = 0; ia < na(is); ia++ )
    {
      D3vector r(rt[is][3*ia+0],rt[is][3*ia+1],rt[is][3*ia+2]);
      r -= rc;
      D3vector v(vt[is][3*ia+0],vt[is][3*ia+1],vt[is][3*ia+2]);
      v -= vc;
      L += mass * ( r ^ v );
    }
  }

  // compute inertia tensor a and vector omega

  // check for special case of all atoms aligned, for which the
  // inertia tensor has a zero eigenvalue

  // check if all atoms are aligned
  // collect all positions in a single vector of D3vectors
  vector<D3vector> rv(size());
  int iat = 0;
  for ( int is = 0; is < vt.size(); is++ )
    for ( int ia = 0; ia < na(is); ia++ )
      rv[iat++] = D3vector(rt[is][3*ia+0],rt[is][3*ia+1],rt[is][3*ia+2]);

  // normalized direction e = (rv[1]-rv[0])
  D3vector e = normalized(rv[1]-rv[0]);
  bool aligned = true;
  for ( int i = 2; (i < size()) && aligned; i++ )
  {
    D3vector u = normalized(rv[i]-rv[0]);
    aligned &= length(u^e) < 1.e-6;
  }

  D3vector omega;
  // compute the inertia tensor
  // treat separately the case of all atoms aligned
  if ( aligned )
  {
    // inertia tensor reduces to a scalar
    double a = 0.0;
    for ( int is = 0; is < vt.size(); is++ )
    {
      double mass = species_list[is]->mass();
      for ( int ia = 0; ia < na(is); ia++ )
      {
        D3vector r(rt[is][3*ia+0],rt[is][3*ia+1],rt[is][3*ia+2]);
        r -= rc;
        a += mass * norm2(r);
      }
    }
    omega = L / a;
  }
  else
  {
    double a00,a01,a02,a10,a11,a12,a20,a21,a22;
    a00=a01=a02=a10=a11=a12=a20=a21=a22=0.0;

    for ( int is = 0; is < vt.size(); is++ )
    {
      double mass = species_list[is]->mass();
      for ( int ia = 0; ia < na(is); ia++ )
      {
        D3vector r(rt[is][3*ia+0],rt[is][3*ia+1],rt[is][3*ia+2]);
        r -= rc;
        a00 += mass * ( r.y * r.y + r.z * r.z );
        a11 += mass * ( r.x * r.x + r.z * r.z );
        a22 += mass * ( r.x * r.x + r.y * r.y );
        a01 -= mass * r.x * r.y;
        a02 -= mass * r.x * r.z;
        a12 -= mass * r.y * r.z;
      }
    }
    a10 = a01;
    a20 = a02;
    a21 = a12;

    // inverse b of the inertia tensor a
    // determinant
    double det = a00 * ( a11 * a22 - a21 * a12 ) -
                 a01 * ( a10 * a22 - a20 * a12 ) +
                 a02 * ( a10 * a21 - a20 * a11 );
    // the determinant must be positive
    assert(det>1.e-8);
    double b00,b01,b02,b10,b11,b12,b20,b21,b22;
    b00 = ( -a12*a21 + a11*a22 ) / det;
    b10 = (  a12*a20 - a10*a22 ) / det;
    b20 = ( -a11*a20 + a10*a21 ) / det;

    b01 = (  a02*a21 - a01*a22 ) / det;
    b11 = ( -a02*a20 + a00*a22 ) / det;
    b21 = (  a01*a20 - a00*a21 ) / det;

    b02 = ( -a02*a11 + a01*a12 ) / det;
    b12 = (  a02*a10 - a00*a12 ) / det;
    b22 = ( -a01*a10 + a00*a11 ) / det;

    // omega = inverse(a) * L
    omega.x = b00 * L.x + b01 * L.y + b02 * L.z;
    omega.y = b10 * L.x + b11 * L.y + b12 * L.z;
    omega.z = b20 * L.x + b21 * L.y + b22 * L.z;
  }

  // correct velocities: v = v - omega ^ r
  for ( int is = 0; is < vt.size(); is++ )
  {
    for ( int ia = 0; ia < na(is); ia++ )
    {
      D3vector r(rt[is][3*ia+0],rt[is][3*ia+1],rt[is][3*ia+2]);
      r -= rc;
      D3vector dv = omega ^ r;
      vt[is][3*ia+0] -= dv.x;
      vt[is][3*ia+1] -= dv.y;
      vt[is][3*ia+2] -= dv.z;
    }
  }
  set_velocities(vt);
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
D3tensor AtomSet::quadrupole(void) const
{
  D3tensor sum;
  for ( int is = 0; is < atom_list.size(); is++ )
  {
    double charge = species_list[is]->zval();
    for ( int ia = 0; ia < atom_list[is].size(); ia++ )
    {
      D3vector p = atom_list[is][ia]->position();
      for ( int idir = 0; idir < 3; idir++ )
        for ( int jdir = 0; jdir < 3; jdir++ )
          sum[idir*3+jdir] += charge * p[idir] * p[jdir];
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::set_cell(const D3vector& a, const D3vector& b, const D3vector& c)
{
  cell_.set(a,b,c);
  sync_cell();
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::sync_cell(void)
{
  sync_cell(cell_);
}

////////////////////////////////////////////////////////////////////////////////
void AtomSet::sync_cell(UnitCell& cell)
{
  double buf[9];
  for ( int i = 0; i < 9; i++ )
    buf[i] = cell.amat(i);
  MPI_Bcast(buf,9,MPI_DOUBLE,0,MPIdata::comm());

  D3vector a0(buf[0],buf[1],buf[2]);
  D3vector a1(buf[3],buf[4],buf[5]);
  D3vector a2(buf[6],buf[7],buf[8]);
  cell.set(a0,a1,a2);
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator << ( ostream &os, const AtomSet &as )
{
  if ( MPIdata::onpe0() )
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
