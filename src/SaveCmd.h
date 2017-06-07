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
// SaveCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAVECMD_H
#define SAVECMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SaveCmd : public Cmd
{
  public:

  Sample *s;

  SaveCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "save"; }
  const char *help_msg(void) const
  {
    return
    "\n save\n\n"
    " syntax: save [-serial] [-text] [-atomsonly] [-no_wfv] filename \n\n"
    "   The save command saves the sample to the file filename.\n\n"
    "   When using the -serial option, I/O is performed from the \n"
    "   head node only. If the -text option is used, wavefunctions\n"
    "   are saved in formatted form instead of base64 encoding. The\n"
    "   -atomsonly option is used to save the atomset only without the\n"
    "   wavefunctions. If the -no_wfv option is used, wavefunction\n"
    "   velocities are not saved.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
