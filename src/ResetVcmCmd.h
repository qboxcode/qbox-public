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
// ResetVcmCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RESETVCMCMD_H
#define RESETVCMCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ResetVcmCmd : public Cmd
{
  public:

  Sample *s;

  ResetVcmCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "reset_vcm"; }
  const char *help_msg(void) const
  {
    return
    "\n reset_vcm\n\n"
    " syntax: reset_vcm \n\n"
    "   The reset_vcm command subtracts the velocity of the center\n"
    "   of mass from the velocity of each atom.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.reset_vcm();
    return 0;
  }
};
#endif
