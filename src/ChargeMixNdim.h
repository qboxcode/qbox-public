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
// ChargeMixNdim.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CHARGEMIXNDIM_H
#define CHARGEMIXNDIM_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class ChargeMixNdim : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "charge_mix_ndim"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("charge_mix_ndim takes one value");

    int v = atoi(argv[1]);
    if ( v < 0 )
      throw invalid_argument("charge_mix_ndim must be non-negative");

    s->ctrl.charge_mix_ndim = v;
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.charge_mix_ndim;
     return st.str();
  }

  ChargeMixNdim(Sample *sample) : s(sample) { s->ctrl.charge_mix_ndim = 3; };
};
#endif
