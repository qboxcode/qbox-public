////////////////////////////////////////////////////////////////////////////////
//
// WfDiag.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WfDiag.h,v 1.3 2007-10-19 16:24:05 fgygi Exp $

#ifndef WFDIAG_H
#define WFDIAG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class WfDiag : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "wf_diag"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " wf_diag takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "T" || v == "F" || v == "EIGVAL" ) )
    {
      if ( ui->onpe0() )
        cout << " wf_diag must be T, F or EIGVAL" << endl;
      return 1;
    }

    s->ctrl.wf_diag = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_diag;
     return st.str();
  }

  WfDiag(Sample *sample) : s(sample) { s->ctrl.wf_diag = "F"; };
};
#endif
