////////////////////////////////////////////////////////////////////////////////
//
// Cell.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Cell.h,v 1.6 2008-02-03 22:53:54 fgygi Exp $

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

  char *name ( void ) const { return "cell"; };

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

    s->atoms.unit_cell.set(a,b,c);

    if ( ui->onpe0() )
    {
      cout << s->atoms.unit_cell;
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
