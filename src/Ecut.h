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
// Ecut.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ECUT_H
#define ECUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecut : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "ecut"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " ecut takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " ecut must be non-negative" << endl;
      return 1;
    }

    if ( s->wf.ecut() == 0.5 * v )
      return 0;

    s->wf.resize(0.5*v);
    if ( s->wfv != 0 )
    {
      s->wfv->resize(0.5*v);
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
     st << setw(10) << 2 * s->wf.ecut();
     return st.str();
  }

  Ecut(Sample *sample) : s(sample) {};
};
#endif
