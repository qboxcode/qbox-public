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
// Nrowmax.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NROWMAX_H
#define NROWMAX_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Nrowmax : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "nrowmax"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " nrowmax takes only one value" << endl;
      return 1;
    }

    int v = atoi(argv[1]);
    if ( v <= 0 )
    {
      if ( ui->onpe0() )
        cout << " nrowmax must be positive" << endl;
      return 1;
    }

    s->wf.set_nrowmax(v);
    s->wf.update_occ(0.0);
    if ( s->wfv != 0 )
    {
      s->wfv->set_nrowmax(v);
      s->wfv->clear();
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nrowmax();
     return st.str();
  }

  Nrowmax(Sample *sample) : s(sample) {};
};
#endif
