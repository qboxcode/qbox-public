////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// AlphaPBE0.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALPHAPBE0_H
#define ALPHAPBE0_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class AlphaPBE0 : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "alpha_PBE0"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("alpha_PBE0 takes one value");

    double v = atof(argv[1]);
    if ( v < 0.0 )
      throw invalid_argument("alpha_PBE0 must be non-negative");

    s->ctrl.alpha_PBE0 = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.alpha_PBE0;
     return st.str();
  }

  AlphaPBE0(Sample *sample) : s(sample)
  {
    s->ctrl.alpha_PBE0 = 0.25;
  }
};
#endif
