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
// Nspin.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NSPIN_H
#define NSPIN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class Nspin : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "nspin"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("nspin takes one value");

    int v = atoi(argv[1]);
    if ( v != 1 && v!=2 )
      throw invalid_argument("nspin must be 1 or 2");

    if ( s->wf.nspin() == v )
      return 0;

    s->wf.set_nspin(v);

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nspin();
     return st.str();
  }

  Nspin(Sample *sample) : s(sample) {};
};
#endif
