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
// PlotCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PLOTCMD_H
#define PLOTCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class PlotCmd : public Cmd
{
  public:

  Sample *s;

  PlotCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "plot"; }
  const char *help_msg(void) const
  {
    return
    "\n plot\n\n"
    " syntax: plot filename\n"
    "         plot -density [-spin {1|2}] filename\n"
    "         plot -vlocal  [-spin {1|2}] filename\n"
    "         plot -wf <n> [-spin {1|2}] filename\n"
    "         plot -wfs <nmin> <nmax> [-spin {1|2}] filename\n\n"
    "   The plot command creates a plot file in xyz or cube format.\n\n"
    "   The default format is xyz, used for plotting atoms only.\n"
    "   When using the -density option, the charge density is written.\n"
    "   When using the -vlocal option, the local potential is written.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
