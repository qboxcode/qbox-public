////////////////////////////////////////////////////////////////////////////////
//
// AtomsDyn.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AtomsDyn.h,v 1.3 2007-10-19 16:24:04 fgygi Exp $

#ifndef ATOMSDYN_H
#define ATOMSDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class AtomsDyn : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "atoms_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " atoms_dyn takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "LOCKED" || v == "SD" || v == "SDA" || v == "MD" ) )
    {
      if ( ui->onpe0() )
        cout << " atoms_dyn must be LOCKED or SD or SDA or MD" << endl;
      return 1;
    }

    s->ctrl.atoms_dyn = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.atoms_dyn;
     return st.str();
  }

  AtomsDyn(Sample *sample) : s(sample) { s->ctrl.atoms_dyn = "LOCKED"; };
};
#endif
