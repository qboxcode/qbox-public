////////////////////////////////////////////////////////////////////////////////
//
// Xc.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Xc.h,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#ifndef XC_H
#define XC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Xc : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "xc"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " xc takes only one value" << endl;
      return 1;
    }
    
    string v = argv[1];
    if ( !( v == "LDA" || v == "PBE" || v == "BLYP" ) )
    {
      if ( ui->onpe0() )
        cout << " xc must be LDA, PBE or BLYP" << endl;
      return 1;
    }

    s->ctrl.xc= v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.xc;
     return st.str();
  }

  Xc(Sample *sample) : s(sample) { s->ctrl.xc = "LDA"; };
};
#endif
