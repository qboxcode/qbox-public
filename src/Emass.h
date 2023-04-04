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
// Emass.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EMASS_H
#define EMASS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class Emass : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "emass"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("emass takes one value");

    double v = atof(argv[1]);
    if ( v < 0.0 )
      throw invalid_argument("emass must be non-negative");

    s->ctrl.emass = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.emass;
     return st.str();
  }

  Emass(Sample *sample) : s(sample) { s->ctrl.emass = 0.0; }
};
#endif
