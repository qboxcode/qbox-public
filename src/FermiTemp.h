////////////////////////////////////////////////////////////////////////////////
//
// FermiTemp.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: FermiTemp.h,v 1.1 2004-10-04 19:55:39 fgygi Exp $

#ifndef FERMITEMP_H
#define FERMITEMP_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class FermiTemp : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "fermi_temp"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " fermi_temp takes only one value" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " fermi_temp must be non-negative" << endl;
      return 1;
    }

    s->ctrl.fermi_temp = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.fermi_temp;
     return st.str();
  }

  FermiTemp(Sample *sample) : s(sample) { s->ctrl.fermi_temp = 0.0; }
};
#endif
