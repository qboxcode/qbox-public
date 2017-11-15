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
// ListSpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LISTSPECIESCMD_H
#define LISTSPECIESCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class ListSpeciesCmd : public Cmd
{
  public:

  Sample *s;

  ListSpeciesCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "list_species"; }
  const char *help_msg(void) const
  {
    return
    "\n list_species\n\n"
    " syntax: list_species\n\n"
    "   The list_species command prints a list of all defined species.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 1 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: list_species" << endl;
      }
      return 1;
    }
    s->atoms.listSpecies();
    return 0;
  }
};
#endif
