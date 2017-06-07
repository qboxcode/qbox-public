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
// FoldInWsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FOLDINWSCMD_H
#define FOLDINWSCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class FoldInWsCmd : public Cmd
{
  public:

  Sample *s;

  FoldInWsCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "fold_in_ws"; }
  const char *help_msg(void) const
  {
    return
    "\n fold_in_ws\n\n"
    " syntax: fold_in_ws \n\n"
    "   The fold_in_ws command folds all atomic positions back in\n"
    "   the Wigner-Seitz cell of the current unit cell.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.fold_in_ws();
    return 0;
  }
};
#endif
