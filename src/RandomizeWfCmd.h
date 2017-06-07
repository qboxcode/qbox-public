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
// RandomizeWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RANDOMIZEWFCMD_H
#define RANDOMIZEWFCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeWfCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeWfCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "randomize_wf"; }
  const char *help_msg(void) const
  {
    return
    "\n randomize_wf\n\n"
    " syntax: randomize_wf [amplitude]\n\n"
    "   The randomize_wf command adds random amplitudes to\n"
    "   the wavefunction Fourier coefficients\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc > 2 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: randomize_wf [amplitude]" << endl;
      }
      return 1;
    }
    double amp = 0.02;
    if ( argc == 2 )
      amp = atof(argv[1]);
    s->wf.randomize(amp);
    return 0;
  }
};
#endif
