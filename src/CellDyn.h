////////////////////////////////////////////////////////////////////////////////
//
// CellDyn.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CellDyn.h,v 1.2 2004-03-11 21:52:31 fgygi Exp $

#ifndef CELLDYN_H
#define CELLDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class CellDyn : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "cell_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " cell_dyn takes only one value" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "LOCKED" || v == "SD" ) )
    {
      if ( ui->onpe0() )
        cout << " cell_dyn must be in [LOCKED,SD]" << endl;
      return 1;
    }

    s->ctrl.cell_dyn = v;
        
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.cell_dyn;
     return st.str();
  }

  CellDyn(Sample *sample) : s(sample) { s->ctrl.cell_dyn = "LOCKED"; };
};
#endif
