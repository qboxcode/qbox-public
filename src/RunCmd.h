////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.h,v 1.4 2007-10-19 16:24:04 fgygi Exp $

#ifndef RUNCMD_H
#define RUNCMD_H

#include <iostream>
#include "UserInterface.h"

class Sample;
class RunCmd : public Cmd
{
  private:

  public:

  Sample *s;

  RunCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "run"; }
  char *help_msg(void) const
  {
    return
    "\n run\n\n"
    " syntax: run n [nite]\n\n"
    "   The run command runs n steps of simulation. Each step\n"
    "   consists of one (optionally nite) electronic steps\n\n";
  }

  int action(int argc, char **argv);

};
#endif
