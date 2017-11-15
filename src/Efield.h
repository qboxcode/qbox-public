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
// Efield.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EFIELD_H
#define EFIELD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include"D3vector.h"

#include "Sample.h"

class Efield : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "e_field"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      if ( ui->onpe0() )
      cout << " e_field takes 3 values" << endl;
      return 1;
    }

    double v0 = atof(argv[1]);
    double v1 = atof(argv[2]);
    double v2 = atof(argv[3]);

    s->ctrl.e_field[0] = v0;
    s->ctrl.e_field[1] = v1;
    s->ctrl.e_field[2] = v2;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.e_field[0] << " "
        << s->ctrl.e_field[1] << " "
        << s->ctrl.e_field[2] << " ";
     return st.str();
  }

  Efield(Sample *sample) : s(sample)
  {
    s->ctrl.e_field[0] = 0.0;
    s->ctrl.e_field[1] = 0.0;
    s->ctrl.e_field[2] = 0.0;
  }
};
#endif
