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
// ExtStress.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EXTSTRESS_H
#define EXTSTRESS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class ExtStress : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "ext_stress"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 7 )
    {
      if ( ui->onpe0() )
      cout << " ext_stress must be specified as s_xx,s_yy,s_zz,s_xy,s_yz,s_xz"
           << endl;
      return 1;
    }

    for ( int i = 0; i < 6; i++ )
      s->ctrl.ext_stress[i] = atof(argv[i+1]);

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     for ( int i = 0; i < 6; i++ )
       st << setw(10) << s->ctrl.ext_stress[i];
     return st.str();
  }

  ExtStress(Sample *sample) : s(sample)
  {
    for ( int i = 0; i < 6; i++ )
      s->ctrl.ext_stress[i] = 0.0;
  };
};
#endif
