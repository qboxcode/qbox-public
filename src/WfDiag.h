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
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class WfDiag : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "wf_diag"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("wf_diag takes one value");

    string v = argv[1];
    if ( !( v == "T" || v == "F" || v == "EIGVAL" ||
            v == "MLWF" || v == "MLWFC" ) )
      throw invalid_argument("wf_diag must be T, F, EIGVAL, MLWF, or MLWFC");

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
