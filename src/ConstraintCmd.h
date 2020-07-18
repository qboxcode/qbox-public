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
// ConstraintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CONSTRAINTCMD_H
#define CONSTRAINTCMD_H

#include "UserInterface.h"
#include "Sample.h"

class ConstraintCmd : public Cmd
{
  public:

  Sample *s;

  ConstraintCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "constraint"; }
  const char *help_msg(void) const
  {
    return
    "\n constraint\n\n"
    " syntax:\n\n"
    "   constraint define position name atom\n"
    "   constraint define distance name atom1 atom2 distance [velocity]\n"
    "   constraint define angle name atom1 atom2 atom3 angle [velocity]\n"
    "   constraint define torsion name atom1 atom2 atom3 atom4 angle [velocity]\n"
    "   constraint set name value [velocity]\n"
    "   constraint delete name\n"
    "   constraint list\n"
    "   constraint enforce\n\n"
    "   Constraints are enforced at each MD step if ions are allowed to move.\n"
    "   If a distance or angle is replaced by '*', the current value is used.\n"
    "   Velocity parameters are optional.\n\n";
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
      return s->constraints.define_constraint(s->atoms,argc,argv);
    }
    else if ( subcmd == "set" )
    {
      return s->constraints.set_constraint(argc,argv);
    }
    else if ( subcmd == "delete" )
    {
      return s->constraints.delete_constraint(argc,argv);
    }
    else if ( subcmd == "enforce" )
    {
      s->constraints.enforce(s->atoms);
      // reset velocities to zero to avoid jump in temperature
      s->atoms.reset_velocities();
      // Note: should catch exception here if constraints could not be enforced
    }
    else if ( subcmd == "list" )
    {
      if ( ui->onpe0() )
        s->constraints.list_constraints(cout);
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
