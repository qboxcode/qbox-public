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
// Thermostat.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Thermostat : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "thermostat"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " thermostat takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "SCALING" || v == "ANDERSEN" || v == "LOWE" ||
         v == "BDP" || v == "OFF" ) )
    {
      if ( ui->onpe0() )
        cout << " thermostat must be SCALING or ANDERSEN or LOWE or BDP or OFF"
             << endl;
      return 1;
    }

    s->ctrl.thermostat = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.thermostat;
     return st.str();
  }

  Thermostat(Sample *sample) : s(sample) { s->ctrl.thermostat = "OFF"; };
};
#endif
