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
// BlHF.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BLHF_H
#define BLHF_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class BlHF : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "blHF"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      if ( ui->onpe0() )
      cout << " blHF takes 3 integer values" << endl;
      return 1;
    }

    int v0 = atoi(argv[1]);
    int v1 = atoi(argv[2]);
    int v2 = atoi(argv[3]);
    if ( v0 < 0 || v1 < 0 || v2 < 0 || v0 > 5 || v1 > 5 || v2 > 5 )
    {
      if ( ui->onpe0() )
        cout << " blHF values must be in [0,5]" << endl;
      return 1;
    }

    s->ctrl.blHF[0] = v0;
    s->ctrl.blHF[1] = v1;
    s->ctrl.blHF[2] = v2;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.blHF[0] << " "
        << s->ctrl.blHF[1] << " "
        << s->ctrl.blHF[2] << " ";
     return st.str();
  }

  BlHF(Sample *sample) : s(sample)
  {
    s->ctrl.blHF[0] = 1;
    s->ctrl.blHF[1] = 1;
    s->ctrl.blHF[2] = 1;
  }
};
#endif
