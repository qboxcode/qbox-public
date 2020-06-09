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
// GlobalExtForce.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "GlobalExtForce.h"
#include "AtomSet.h"
#include "Atom.h"
#include "Species.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void GlobalExtForce::setup(const AtomSet& atoms)
{
  // All atoms are affected, no setup needed
  na_ = atoms.size();
}

////////////////////////////////////////////////////////////////////////////////
double GlobalExtForce::energy(const vector<vector<double> > &r,
                              vector<vector<double> > &f)
{
  double sum = 0.0;
  for ( int is = 0; is < f.size(); is++ )
  {
    const int nais = f[is].size() / 3;
    for ( int ia = 0; ia < nais; ia++ )
    {
      f[is][3*ia+0] += force_.x / na_;
      f[is][3*ia+1] += force_.y / na_;
      f[is][3*ia+2] += force_.z / na_;
      sum -= force_.x * r[is][3*ia+0] / na_ +
             force_.y * r[is][3*ia+1] / na_ +
             force_.z * r[is][3*ia+2] / na_;
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
ostream& GlobalExtForce::print( ostream &os )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <extforce name=\"" << name();
  os << "\" type=\"" << type();
  os << "\" atoms=\"_all_\"\n";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << "            value=\"" << setprecision(6) << force_.x << " "
     << setprecision(6) << force_.y << " "
     << setprecision(6) << force_.z << "\"/>";
  return os;
}
