////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.h,v 1.1 2003-01-25 01:23:31 fgygi Exp $

#ifndef RUNCMD_H
#define RUNCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"

class RunCmd : public Cmd
{
  private:

  int niter, nfi;

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
