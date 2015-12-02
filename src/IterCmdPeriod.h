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
// IterCmdPeriod.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ITERCMDPERIOD_H
#define ITERCMDPERIOD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class IterCmdPeriod : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "iter_cmd_period"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " iter_cmd_period must be a positive integer" << endl;
      return 1;
    }

    int v = atoi(argv[1]);
    if ( v <= 0 )
    {
      if ( ui->onpe0() )
        cout << " iter_cmd_period must be a positive integer" << endl;
      return 1;
    }
    s->ctrl.iter_cmd_period = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.iter_cmd_period;
     return st.str();
  }

  IterCmdPeriod(Sample *sample) : s(sample) { s->ctrl.iter_cmd_period = 1; };
};
#endif
