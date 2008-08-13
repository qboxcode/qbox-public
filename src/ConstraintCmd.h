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
// $Id: ConstraintCmd.h,v 1.4 2008-08-13 06:39:42 fgygi Exp $

#ifndef CONSTRAINTCMD_H
#define CONSTRAINTCMD_H

#include "UserInterface.h"
#include "Sample.h"

class ConstraintCmd : public Cmd
{
  public:

  Sample *s;

  ConstraintCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "constraint"; }
  char *help_msg(void) const
  {
    return
    "\n constraint\n\n"
    " syntax:\n\n"
    "   constraint define distance name atom1 atom2 distance [velocity]\n"
    "   constraint define angle name atom1 atom2 atom3 angle [velocity]\n"
    "   constraint define torsion name atom1 atom2 atom3 atom4 angle [velocity]\n"
    "   constraint set name value [velocity]\n"
    "   constraint delete name\n"
    "   constraint list\n"
    "   constraint enforce\n\n"
    "   Constraints are enforced at each MD step if ions are allowed to move.\n"
    "   Velocity parameters are optional.\n\n";
    //"   constraint multidistance weight name1 name2 \n"
    //"                           [weight name1 name2] distance [velocity]\n"
  }

  int action(int argc, char **argv)
  {
    const bool onpe0 = s->ctxt_.onpe0();
    if ( argc < 2 )
    {
      if ( onpe0 )
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
      // Note: should catch exception here if constraints could not be enforced
    }
    else if ( subcmd == "list" )
    {
      if ( onpe0 )
        s->constraints.list_constraints(cout);
    }
    else
    {
      if ( onpe0 )
        cout << help_msg();
    }

    return 0;
  }
};
#endif
