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
// Debug.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Debug.h,v 1.4 2008-08-13 06:39:42 fgygi Exp $

#ifndef DEBUG_H
#define DEBUG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Debug : public Var
{
  Sample *s;

  public:

  char *name ( void ) const { return "debug"; };

  int set ( int argc, char **argv )
  {
    string v;
    for ( int iarg = 1; iarg < argc; iarg++ )
    {
      string vt = argv[iarg];
      v += " " + vt;
    }

    s->ctrl.debug = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.debug;
     return st.str();
  }

  Debug(Sample *sample) : s(sample) { s->ctrl.debug = "OFF"; }
};
#endif
