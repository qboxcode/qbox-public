////////////////////////////////////////////////////////////////////////////////
//
// ListConstraintsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ListConstraintsCmd.h,v 1.1 2005-06-27 22:34:46 fgygi Exp $

#ifndef LISTCONSTRAINTSCMD_H
#define LISTCONSTRAINTSCMD_H

#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ListConstraintsCmd : public Cmd
{
  public:

  Sample *s;

  ListConstraintsCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "list_constraints"; }
  char *help_msg(void) const
  {
    return 
    "\n list_constraints\n\n"
    " syntax: list_constraints\n\n"
    "   The list_constraints command prints the list"
    " of active constraints.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( s->ctxt_.onpe0() ) s->constraints.list_constraints(cout);
    return 0;
  }
};
#endif
