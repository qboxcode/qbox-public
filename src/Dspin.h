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
#include<stdlib.h>

#include "Sample.h"

class Dspin : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "delta_spin"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
        cout << " delta_spin takes only one value" << endl;
      return 1;
    }

    int v = atoi(argv[1]);
    if ( v < 0 )
    {
      if ( ui->onpe0() )
        cout << " delta_spin must >= 0" << endl;
      return 1;
    }

    if ( s->wf.nspin() < 2 )
    {
      if ( ui->onpe0() )
        cout << " Cannot set delta_spin: nspin == 1" << endl;
      return 1;
    }

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
