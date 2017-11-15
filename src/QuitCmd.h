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
// QuitCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef QUITCMD_H
#define QUITCMD_H

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>

#include "UserInterface.h"
#include "Sample.h"

class QuitCmd : public Cmd
{
  public:

  Sample* s;

  QuitCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "quit"; }
  const char *help_msg(void) const
  {
    return
    "\n quit\n\n"
    " syntax: quit\n\n"
    "   The quit command exits without saving any data.\n\n";
  }

  int action(int argc, char **argv)
  {
    ui->terminate();
    return 0;
  }
};
#endif
