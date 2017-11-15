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
// RunCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

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

  const char *name(void) const { return "run"; }
  const char *help_msg(void) const
  {
    return
    "\n run\n\n"
    " syntax: run [-atomic_density] n [nitscf [nite]]\n\n"
    "   The run command runs n steps of simulation. Each step\n"
    "   consists of one or more (nitscf) scf steps, each consisting\n"
    "   of one or more (nite) electronic steps.\n"
    "   If the -atomic_density option is used, the initial charge\n"
    "   density is a sum of atomic charge densities.\n\n";
  }

  int action(int argc, char **argv);

};
#endif
