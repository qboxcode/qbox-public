////////////////////////////////////////////////////////////////////////////////
//
// Dt.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Dt.h,v 1.1 2003-06-11 22:10:11 fgygi Exp $

#ifndef DT_H
#define DT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Dt : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "dt"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " dt takes only one value" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " dt must be non-negative" << endl;
      return 1;
    }

    s->ctrl.dt = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.dt;
     return st.str();
  }

  Dt(Sample *sample) : s(sample) { s->ctrl.dt = 3.0; }
};
#endif
