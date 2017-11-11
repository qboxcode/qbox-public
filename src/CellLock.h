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
// CellLock.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CELLLOCK_H
#define CELLLOCK_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class CellLock : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "cell_lock"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " cell_lock takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "OFF" || v == "A" || v == "B" || v == "C" ||
            v == "AB" || v == "AC" || v == "BC" || v == "ABC" ||
            v == "S"  || v == "AS" || v == "BS" || v == "CS" ||
            v == "ABS" || v == "ACS" || v == "BCS" || v == "R") )
    {
      if ( ui->onpe0() )
        cout << " cell_lock must be in "
             << "[OFF,A,B,C,AB,AC,BC,S,AS,BS,CS,ABS,ACS,BCS,R]" << endl;
      return 1;
    }

    s->ctrl.cell_lock = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_lock;
     return st.str();
  }

  CellLock(Sample *sample) : s(sample) { s->ctrl.cell_lock = "OFF"; }
};
#endif
