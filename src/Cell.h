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
// Cell.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CELL_H
#define CELL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Cell : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "cell"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 10 )
    {
      if ( ui->onpe0() )
      cout << " cell must be specified with 3 vectors (9 values)" << endl;
      return 1;
    }

    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);

    if ( cell.volume() < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " cell volume must be positive" << endl;
      return 1;
    }

    s->wf.resize(cell,s->wf.refcell(),s->wf.ecut());
    if ( s->wfv != 0 )
    {
      s->wfv->resize(cell,s->wf.refcell(),s->wf.ecut());
      s->wfv->clear();
    }

    s->atoms.set_cell(a,b,c);

    if ( ui->onpe0() )
    {
      cout << s->atoms.cell();
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.cell();
     return st.str();
  }

  Cell(Sample *sample) : s(sample) {};
};
#endif
