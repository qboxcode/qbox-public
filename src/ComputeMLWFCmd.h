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
// ComputeMLWFCmd.h
//
////////////////////////////////////////////////////////////////////////////////

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

  public:
  Sample *s;

  ComputeMLWFCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "compute_mlwf"; }
  const char *help_msg(void) const
  {
    return
    "\n compute_mlwf\n\n"
    " syntax: compute_mlwf\n\n"
    "   The compute_mlwf command computes maximally localized \n"
    " Wannier functions.\n\n";
  }

  int action(int argc, char **argv);

  ComputeMLWFCmd();
};
#endif
