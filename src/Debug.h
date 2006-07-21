////////////////////////////////////////////////////////////////////////////////
//
// Debug.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Debug.h,v 1.2 2006-07-21 17:48:56 fgygi Exp $

#ifndef DEBUG_H
#define DEBUG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Debug : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "debug"; };

  int set ( int argc, char **argv )
  {
    string v;
    for ( int iarg = 1; iarg < argc; iarg++ )
    {
      string vt = argv[iarg];
      v += " " + vt;
    }
    
    s->ctrl.debug = v;
    
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.debug;
     return st.str();
  }

  Debug(Sample *sample) : s(sample) { s->ctrl.debug = "OFF"; }
};
#endif
