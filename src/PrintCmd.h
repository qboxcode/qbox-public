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
// PrintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PRINTCMD_H
#define PRINTCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"

class PrintCmd : public Cmd
{
  public:

  Sample *s;

  PrintCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "print"; }
  const char *help_msg(void) const
  {
    return
    "\n print\n\n"
    " syntax: print variable\n\n"
    "   The print command prints the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
        cout << " use: print variable" << endl;
      return 1;
    }

    Var *varptr = ui->findVar(argv[1]);
    if ( varptr )
    {
      string s = varptr->print();
      if ( ui->onpe0() )
        cout << s << endl;
    }
    else
    {
      // variable is not in the variable list
      if ( ui->onpe0() )
      {
        cout << " no such variable: " << argv[1] << endl;
      }
      return 1;
    }
    return 0;
  }
};
#endif
