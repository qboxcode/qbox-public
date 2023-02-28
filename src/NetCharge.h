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
// NetCharge.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NETCHARGE_H
#define NETCHARGE_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<cstdlib>
#include<stdexcept>

#include "Sample.h"

class NetCharge : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "net_charge"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
      throw invalid_argument("net_charge takes one value");

    const int v = atoi(argv[1]);

    // compute the current netcharge
    // Definition: wf_nel = atoms_nel - netcharge
    // Note: the sign of netcharge is negative if extra electrons are present
    const int netcharge_before = s->atoms.nel() - s->wf.nel();
    if ( v == netcharge_before )
      return 0;

    // set new netcharge to v
    if ( s->atoms.nel() - v < 0 )
      throw invalid_argument("net_charge: cannot remove more than " +
                              to_string(s->atoms.nel()) + " electrons");

    s->wf.set_nel(s->atoms.nel() - v);
    s->wf.update_occ(0.0);
    if ( s->wfv != 0 )
    {
      s->wfv->set_nel(s->atoms.nel() - v);
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->atoms.nel() - s->wf.nel();
     return st.str();
  }

  NetCharge(Sample *sample) : s(sample) {};
};
#endif
