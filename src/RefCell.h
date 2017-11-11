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
// RefCell.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef REFCELL_H
#define REFCELL_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class RefCell : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "ref_cell"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 10 )
    {
      if ( ui->onpe0() )
      cout << " ref_cell must be specified with 3 vectors (9 values)" << endl;
      return 1;
    }

    D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell ref_cell(a0,a1,a2);

    if ( ref_cell.volume() < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " ref_cell volume must be positive" << endl;
      return 1;
    }

    s->wf.resize(s->wf.cell(), ref_cell,s->wf.ecut());
    if ( s->wfv != 0 )
    {
      s->wfv->resize(s->wf.cell(),s->wf.refcell(),s->wf.ecut());
      s->wfv->clear();
    }

    if ( ui->onpe0() )
    {
      cout << "  <reference_unit_cell>\n"
           << s->wf.refcell()
           << "  </reference_unit_cell>" << endl;
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.refcell();
     return st.str();
  }

  RefCell(Sample *sample) : s(sample) {}
};
#endif
