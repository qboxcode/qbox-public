////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
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
// ExtForceCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EXTFORCECMD_H
#define EXTFORCECMD_H

#include "UserInterface.h"
#include "Sample.h"

class ExtForceCmd : public Cmd
{
  public:

  Sample *s;

  ExtForceCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "extforce"; }
  const char *help_msg(void) const
  {
    return
    "\n extforce\n\n"
    " syntax:\n\n"
    "   extforce define atomic name atom fx fy fz\n"
    "   extforce define pair name atom1 atom2 force\n"
    "   extforce define global name fx fy fz\n"
    "   extforce set name fx fy fz\n"
    "   extforce set name f\n"
    "   extforce delete name\n"
    "   extforce list\n"
    "   External forces are added to ionic forces at each MD step.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc < 2 )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    string subcmd(argv[1]);
    if ( subcmd == "define" )
    {
      return s->extforces.define_extforce(s->atoms,argc,argv);
    }
    else if ( subcmd == "set" )
    {
      return s->extforces.set_extforce(argc,argv);
    }
    else if ( subcmd == "delete" )
    {
      return s->extforces.delete_extforce(argc,argv);
    }
    else if ( subcmd == "list" )
    {
      if ( ui->onpe0() )
        s->extforces.list_extforces(cout);
    }
    else
    {
      if ( ui->onpe0() )
        cout << help_msg();
    }

    return 0;
  }
};
#endif
