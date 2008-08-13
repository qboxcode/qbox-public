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
// Xc.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Xc.h,v 1.4 2008-08-13 06:39:43 fgygi Exp $

#ifndef XC_H
#define XC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Xc : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "xc"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " xc takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "LDA" || v == "PBE" || v == "BLYP" ) )
    {
      if ( ui->onpe0() )
        cout << " xc must be LDA, PBE or BLYP" << endl;
      return 1;
    }

    s->ctrl.xc= v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.xc;
     return st.str();
  }

  Xc(Sample *sample) : s(sample) { s->ctrl.xc = "LDA"; };
};
#endif
