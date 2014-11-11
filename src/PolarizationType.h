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
// PolarizationType.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef POLARIZATION_TYPE_H
#define POLARIZATION_TYPE_H

#include<iostream>
#include<iomanip>
#include<sstream>

#include "Sample.h"

class PolarizationType : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "polarization_type"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " polarization_type takes only one value" << endl;
      return 1;
    }

    string v = argv[1];

    if ( v == "MLWF" || v == "MLWF_REF" || v == "BERRY" )
      s->ctrl.polarization_type = v;
    else
    {
      if ( ui->onpe0() )
      cout << " polarization_type must be MLWF, MLWF_REF or BERRY" << endl;
      return 1;
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.polarization_type;
     return st.str();
  }

  PolarizationType(Sample *sample) : s(sample)
  {
    s->ctrl.polarization_type = "MLWF";
  }
};
#endif
