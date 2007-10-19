////////////////////////////////////////////////////////////////////////////////
//
// ChargeMixRcut.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeMixRcut.h,v 1.2 2007-10-19 16:24:04 fgygi Exp $

#ifndef CHARGEMIXRCUT_H
#define CHARGEMIXRCUT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ChargeMixRcut : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "charge_mix_rcut"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " charge_mix_rcut takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " charge_mix_rcut must be non-negative" << endl;
      return 1;
    }
    s->ctrl.charge_mix_rcut = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.charge_mix_rcut;
     return st.str();
  }

  ChargeMixRcut(Sample *sample) : s(sample) { s->ctrl.charge_mix_rcut = 10.0; };
};
#endif
