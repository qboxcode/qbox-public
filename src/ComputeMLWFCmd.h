////////////////////////////////////////////////////////////////////////////////
//
// ComputeMLWFCmd.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ComputeMLWFCmd.h,v 1.2 2007-10-19 16:24:04 fgygi Exp $

#ifndef COMPUTEMLWFCMD_H
#define COMPUTEMLWFCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"
#include "MLWFTransform.h"

class ComputeMLWFCmd : public Cmd
{
  private:

  MLWFTransform* mlwft;

  public:
  Sample *s;

  ComputeMLWFCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "compute_mlwf"; }

  char *help_msg(void) const
  {
    return
    "\n compute_mlwf\n\n"
    " syntax: compute_mlwf\n\n"
    "   The compute_mlwf command computes maximally localized \n"
    " Wannier functions.\n\n";
  }

  int action(int argc, char **argv);

  ComputeMLWFCmd();
  ~ComputeMLWFCmd();
};
#endif
