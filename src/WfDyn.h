////////////////////////////////////////////////////////////////////////////////
//
// WfDyn.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WfDyn.h,v 1.6 2007-10-19 16:24:05 fgygi Exp $

#ifndef WFDYN_H
#define WFDYN_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"

class WfDyn : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "wf_dyn"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " wf_dyn takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "LOCKED" || v == "SD" || v == "PSD" ||
            v == "PSDA" || v == "MD" ) )
    {
      if ( ui->onpe0() )
        cout << " wf_dyn must be in [LOCKED,SD,PSD,PSDA,MD]" << endl;
      return 1;
    }

    s->ctrl.wf_dyn = v;

//     if ( v == "MD" )
//     {
//       if ( s->wfv == 0 )
//       {
//         s->wfv = new Wavefunction(s->wf);
//         s->wfv->clear();
//       }
//     }
//     else
//     {
//       delete s->wfv;
//       s->wfv = 0;
//     }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_dyn;
     return st.str();
  }

  WfDyn(Sample *sample) : s(sample) { s->ctrl.wf_dyn = "SD"; };
};
#endif
