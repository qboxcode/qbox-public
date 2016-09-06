////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2015 The Regents of the University of California
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
// IterCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ITERCMD_H
#define ITERCMD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class IterCmd : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "iter_cmd"; };

  int set ( int argc, char **argv )
  {
    if ( !strcmp(argv[1],"NULL") )
    {
      // reset iter_cmd to empty string
      s->ctrl.iter_cmd.clear();
      return 0;
    }

    // include all arguments until \n
    string v;
    for ( int i = 1; i < argc; i++ )
    {
      if ( i > 1 )
        v += " ";
      v += argv[i];
    }

    s->ctrl.iter_cmd = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.iter_cmd;
     return st.str();
  }

  IterCmd(Sample *sample) : s(sample) { s->ctrl.iter_cmd = ""; };
};
#endif
