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
// $Id: SaveCmd.h,v 1.4 2008-08-13 06:39:43 fgygi Exp $

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

  char *name(void) const { return "save"; }

  char *help_msg(void) const
  {
    return
    "\n save\n\n"
    " syntax: save [-serial] [-text] filename \n\n"
    "   The save command saves the sample to the file filename.\n\n"
    "   When using the -serial option, I/O is performed from the \n"
    "   head node only. \n\n";
  }

  int action(int argc, char **argv);
};
#endif
