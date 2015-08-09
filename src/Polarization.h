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
// Polarization.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef POLARIZATION_H
#define POLARIZATION_H

#include<iostream>
#include<iomanip>
#include<sstream>

#include "Sample.h"

class Polarization: public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "polarization"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " polarization takes only one value" << endl;
      return 1;
    }

    string v = argv[1];

    if ( v == "OFF" || v == "MLWF" || v == "MLWF_REF" || v == "MLWF_REF_Q" ||
         v == "BERRY" )
      s->ctrl.polarization = v;
    else
    {
      if ( ui->onpe0() )
      cout <<
      " polarization must be OFF, MLWF, MLWF_REF, MLWF_REF_Q or BERRY" << endl;
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
     st << s->ctrl.polarization;
     return st.str();
  }

  Polarization(Sample *sample) : s(sample)
  {
    s->ctrl.polarization = "OFF";
  }
};
#endif
