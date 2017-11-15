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
// LoadCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef LOADCMD_H
#define LOADCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class LoadCmd : public Cmd
{
  public:

  Sample *s;

  LoadCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "load"; }
  const char *help_msg(void) const
  {
    return
    "\n load\n\n"
    " syntax: load [-serial] filename \n\n"
    "   The load command loads a sample from the file filename.\n"
    "   The -serial option bypasses the parallel load algorithm.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
