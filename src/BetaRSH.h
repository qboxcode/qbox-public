////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// BetaRSH.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BETARSH_H
#define BETARSH_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class BetaRSH : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "beta_RSH"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " beta_RSH takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " beta_RSH must be non-negative" << endl;
      return 1;
    }

    s->ctrl.beta_RSH = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.beta_RSH;
     return st.str();
  }

  BetaRSH(Sample *sample) : s(sample)
  {
    s->ctrl.beta_RSH = 0.25;
  }
};
#endif
