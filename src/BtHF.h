////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// BtHF.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BTHF_H
#define BTHF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class BtHF : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "btHF"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " btHF takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v >= 1.0 || v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " btHF value must be in [0,1)" << endl;
      return 1;
    }

    s->ctrl.btHF = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.btHF;
     return st.str();
  }

  BtHF(Sample *sample) : s(sample)
  {
    s->ctrl.btHF = 0.0;
  }
};
#endif
