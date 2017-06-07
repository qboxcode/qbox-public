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
// ListAtomsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LISTATOMSCMD_H
#define LISTATOMSCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class ListAtomsCmd : public Cmd
{
  public:

  Sample *s;

  ListAtomsCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "list_atoms"; }
  const char *help_msg(void) const
  {
    return
    "\n list_atoms\n\n"
    " syntax: list_atoms\n\n"
    "   The list_atoms command prints a list of all defined atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: list_atoms" << endl;
      }
      return 1;
    }
    s->atoms.listAtoms();
    return 0;
  }
};
#endif
