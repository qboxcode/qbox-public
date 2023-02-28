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
// AtomsDyn.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ATOMSDYN_H
#define ATOMSDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class AtomsDyn : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "atoms_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("atoms_dyn takes one value");

    string v = argv[1];
    if ( !( v == "LOCKED" ||
            v == "SD" ||
            v == "SDA" ||
            v == "CG" ||
            v == "AND" ||
            v == "MD" ||
            v == "BMD" ) )
    throw invalid_argument("atoms_dyn must be "
                           "LOCKED, SD, SDA, CG, AND, MD or BMD]");
    s->ctrl.atoms_dyn = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.atoms_dyn;
     return st.str();
  }

  AtomsDyn(Sample *sample) : s(sample) { s->ctrl.atoms_dyn = "LOCKED"; };
};
#endif
