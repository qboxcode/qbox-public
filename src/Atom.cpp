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
// Atom.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Atom.h"
#include <iomanip>
using namespace std;

Atom::Atom (string newname, string newspecies, D3vector pos, D3vector vel)
{
  name_ = newname;
  species_ = newspecies;
  position_ = pos;
  velocity_ = vel;
}

ostream& operator << ( ostream &os, Atom &a )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <atom name=\"" << a.name() << "\""
     << " species=\"" << a.species() << "\">\n"
     << "    <position> ";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setw(12) << setprecision(8) << a.position().x << " "
     << setw(12) << setprecision(8) << a.position().y << " "
     << setw(12) << setprecision(8) << a.position().z << "  "
     << " </position>\n"
     << "    <velocity> ";
  os.setf(ios::scientific,ios::floatfield);
  os << setw(13) << setprecision(6) << a.velocity().x << " "
     << setw(13) << setprecision(6) << a.velocity().y << " "
     << setw(13) << setprecision(6) << a.velocity().z
     << " </velocity>\n  </atom>" << endl;
  return os;
}
