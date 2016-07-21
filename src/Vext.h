////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2016 The Regents of the University of California
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
// Vext.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VEXT_H
#define VEXT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "ExternalPotential.h"

class Vext : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "vext"; };

  int set ( int argc, char **argv )
  {
    if ( argc > 3 )
    {
      if ( ui->onpe0() )
      cout << " vext takes only one value" << endl;
      return 1;
    }

    if ( !strcmp(argv[1],"NULL") )
    {
      // set vext NULL
      // reset file name to empty string
      delete s->vext;
      s->vext = 0;
      s->ctrl.vext.clear();
    }
    else
    {
      s->ctrl.vext = argv[1];
      s->vext = new ExternalPotential(*s);
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.vext;
     return st.str();
  }

  Vext(Sample *sample) : s(sample) { s->ctrl.vext = ""; };
};
#endif
