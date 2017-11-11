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
// WfDiag.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef WFDIAG_H
#define WFDIAG_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class WfDiag : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "wf_diag"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " wf_diag takes only one value" << endl;
      return 1;
    }

    string v = argv[1];
    if ( !( v == "T" || v == "F" || v == "EIGVAL" ||
            v == "MLWF" || v == "MLWFC" ) )
    {
      if ( ui->onpe0() )
        cout << " wf_diag must be in T, F, EIGVAL, MLWF, MLWFC" << endl;
      return 1;
    }

    s->ctrl.wf_diag = v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.wf_diag;
     return st.str();
  }

  WfDiag(Sample *sample) : s(sample) { s->ctrl.wf_diag = "F"; };
};
#endif
