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
// ResetRotationCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RESETROTATIONCMD_H
#define RESETROTATIONCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class ResetRotationCmd : public Cmd
{
  public:

  Sample *s;

  ResetRotationCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "reset_rotation"; }
  const char *help_msg(void) const
  {
    return
    "\n reset_rotation\n\n"
    " syntax: reset_rotation\n\n"
    "   The reset_rotation command adjusts velocities of the atoms\n"
    "   to cancel the global rotation of the system.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.reset_rotation();
    return 0;
  }
};
#endif
