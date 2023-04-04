////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012 The Regents of the University of California
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
// Dspin.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef DSPIN_H
#define DSPIN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class Dspin : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "delta_spin"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("delta_spin takes one value");

    int v = atoi(argv[1]);
    if ( v < 0 )
      throw invalid_argument("delta_spin must be non-negative");

    if ( s->wf.nspin() < 2 )
      throw invalid_argument("Cannot set delta_spin, nspin==1");

    s->wf.set_deltaspin(v);

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.deltaspin();
     return st.str();
  }

  Dspin(Sample *sample) : s(sample) {};
};
#endif
