////////////////////////////////////////////////////////////////////////////////
//
// Nrowmax.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Nrowmax.h,v 1.1 2003-01-25 01:23:31 fgygi Exp $

#ifndef NROWMAX_H
#define NROWMAX_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Nrowmax : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "nrowmax"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " nrowmax takes only one value" << endl;
      return 1;
    }
    
    int v = atoi(argv[1]);
    if ( v <= 0 )
    {
      if ( ui->onpe0() )
        cout << " nrowmax must be positive" << endl;
      return 1;
    }

    s->wf.set_nrowmax(v);
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->wf.nrowmax();
     return st.str();
  }

  Nrowmax(Sample *sample) : s(sample) {};
};
#endif
