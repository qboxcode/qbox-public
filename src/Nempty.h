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
// Nempty.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEMPTY_H
#define NEMPTY_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Nempty : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "nempty"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " nempty takes only one value" << endl;
      return 1;
    }

    int v = atoi(argv[1]);
    if ( v < 0 )
    {
      if ( ui->onpe0() )
        cout << " nempty must be non-negative" << endl;
      return 1;
    }

    s->wf.set_nempty(v);
    if ( s->wfv != 0 )
    {
      s->wfv->set_nempty(v);
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nempty();
     return st.str();
  }

  Nempty(Sample *sample) : s(sample) {};
};
#endif
