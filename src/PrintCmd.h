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
// $Id: PrintCmd.h,v 1.3 2008-08-13 06:39:43 fgygi Exp $

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

  char *name(void) const { return "print"; }

  char *help_msg(void) const
  {
    return
    "\n print\n\n"
    " syntax: print variable\n\n"
    "   The print command prints the value of an interface variable.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( ui->onpe0() )
    {
      if ( argc != 2 )
      {
        cout << " use: print variable" << endl;
        return 1;
      }

      Var *varptr = ui->findVar(argv[1]);
      if ( varptr )
      {
        cout << varptr->print() << endl;
      }
      else
      {
        // variable is not in the variable list
        cout << " no such variable: " << argv[1] << endl;
        return 1;
      }
    }
    return 0;
  }
};
#endif
