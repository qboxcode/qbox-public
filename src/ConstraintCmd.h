////////////////////////////////////////////////////////////////////////////////
//
// ConstraintCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintCmd.h,v 1.1 2005-06-27 22:34:46 fgygi Exp $

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
    "   constraint distance name1 name2 distance [velocity]\n"
    "   The constraint command defines a new constraint and adds it to\n"
    "   the constraint list. Constraints are enforced at each MD step\n"
    "   if ions are allowed to move. The optional velocity parameter\n"
    "   must be given in atomic units.\n\n";
    //"   constraint multidistance weight name1 name2 \n"
    //"                           [weight name1 name2] distance [velocity]\n"
  }

  int action(int argc, char **argv)
  {
    string subcmd(argv[1]);
    if ( subcmd == "distance" || subcmd == "angle" || subcmd == "torsion" )
      return s->constraints.set_constraint(s->atoms,argc,argv);
    else if ( subcmd == "enforce" )
    {
      s->constraints.enforce(s->atoms);
      // Note: should catch exception here if constraints could not be enforced
    }
    return 0;
  }
};
#endif
