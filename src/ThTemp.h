////////////////////////////////////////////////////////////////////////////////
//
// ThTemp.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ThTemp.h,v 1.1 2003-08-28 00:28:51 fgygi Exp $

#ifndef THTEMP_H
#define THTEMP_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ThTemp : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "th_temp"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " th_temp takes only one value" << endl;
      return 1;
    }
    
    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " th_temp must be non-negative" << endl;
      return 1;
    }

    s->ctrl.th_temp = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.th_temp;
     return st.str();
  }

  ThTemp(Sample *sample) : s(sample) { s->ctrl.th_temp = 0.0; }
};
#endif
