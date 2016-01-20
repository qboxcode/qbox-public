////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2015 The Regents of the University of California
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
// PartialChargeCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PARTIALCHARGECMD_H
#define PARTIALCHARGECMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class PartialChargeCmd : public Cmd
{
  public:

  Sample *s;

  PartialChargeCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "partial_charge"; }
  const char *help_msg(void) const
  {
    return
    "\n partial_charge\n\n"
    " syntax: partial_charge [-spin {1|2}] name radius\n"
    "   The partial_charge command computes the amount of charge\n"
    "   density contained in a sphere centered on an atom.\n"
    "   When using the -spin option, the charge of the given spin\n"
    "   is computed.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
