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
// SetCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SETCMD_H
#define SETCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SetCmd : public Cmd
{
  public:

  Sample *s;

  SetCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "set"; }
  const char *help_msg(void) const
  {
    return
    "\n set\n\n"
    " syntax: set variable value[s]\n\n"
    "   The set command sets the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc < 3 )
    {
      if ( ui->onpe0() )
        cout << " use: set variable value[s]" << endl;
      return 1;
    }

    Var* varptr = ui->findVar(argv[1]);

    if ( varptr )
    {
      varptr->set(argc-1,&argv[1]);
    }
    else
    {
      // variable is not in the variable list
      if ( ui->onpe0() )
        cout << " no such variable: " << argv[1] << endl;
      return 1;
    }
    return 0;
  }
};
#endif
