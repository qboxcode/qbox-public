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

  const char *name ( void ) const { return "debug"; };

  int set ( int argc, char **argv )
  {
    // use: set debug key [val ..]
    // use: set debug key
    if ( argc < 2 )
    {
      if ( ui->onpe0() )
      cout << " use: set debug key val [val ...]" << endl;
      cout << " use: set debug key" << endl;
      return 1;
    }
    string key(argv[1]);
    // if ( ui->onpe0() ) cout << "Debug: key = " << key << endl;
    string val;
    for ( int iarg = 2; iarg < argc; iarg++ )
    {
      string vt = argv[iarg];
      if ( iarg > 2 )
        val += " ";
      val += vt;
    }
    // if ( ui->onpe0() ) cout << "Debug: val = " << val << endl;

    if ( val.empty() )
    {
      if ( ui->onpe0() ) cout << "Debug: reset key " << key << endl;
      s->ctrl.debug.erase(key);
    }
    else
    {
      // add key,value pair to debug map
      s->ctrl.debug[key] = val;
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = " << endl;
     st.setf(ios::right,ios::adjustfield);
     for ( std::map<std::string,std::string>::iterator i =
           s->ctrl.debug.begin(); i != s->ctrl.debug.end(); ++i )
     {
       st << setw(12) << i->first << " " << i->second << endl;
     }
     return st.str();
  }

  Debug(Sample *sample) : s(sample) {}
};
#endif
