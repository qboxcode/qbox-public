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
// SpectrumCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SPECTRUMCMD_H
#define SPECTRUMCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"
#include "MLWFTransform.h"

class SpectrumCmd : public Cmd
{
  private:

  public:
  Sample *s;

  SpectrumCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "spectrum"; }
  const char *help_msg(void) const
  {
    return
    "\n spectrum\n\n"
    " syntax: spectrum\n\n"
    "   The spectrum command computes dipole matrix elements \n"
    " between occupied and empty orbitals.\n\n";
  }

  int action(int argc, char **argv);

  SpectrumCmd();
};
#endif
