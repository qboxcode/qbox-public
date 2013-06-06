////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// RandomizeRCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RANDOMIZERCMD_H
#define RANDOMIZERCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeRCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeRCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "randomize_r"; }
  const char *help_msg(void) const
  {
    return
    "\n randomize_r\n\n"
    " syntax: randomize_r amplitude\n\n"
    "   The randomize_r command adds random displacements to all atoms\n"
    "   using random numbers drawn from a normal distribution\n"
    "   scaled by the amplitude parameter.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: randomize_r amplitude" << endl;
      }
      return 1;
    }
    const double amplitude = atof(argv[1]);
    if ( amplitude < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " randomize_r: amplitude must be non-negative" << endl;
      return 1;
    }

    s->atoms.randomize_positions(amplitude);
    return 0;
  }
};
#endif
