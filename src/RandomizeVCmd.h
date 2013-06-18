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
// RandomizeVCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RANDOMIZEVCMD_H
#define RANDOMIZEVCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class RandomizeVCmd : public Cmd
{
  public:

  Sample *s;

  RandomizeVCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "randomize_v"; }
  const char *help_msg(void) const
  {
    return
    "\n randomize_v\n\n"
    " syntax: randomize_v temp\n\n"
    "   The randomize_v command initializes atomic velocities\n"
    "   with random numbers drawn from a Maxwell-Boltzmann distribution\n"
    "   at temperature temp\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: randomize_v temp" << endl;
      }
      return 1;
    }
    const double temp = atof(argv[1]);
    if ( temp <= 0.0 )
    {
      if ( ui->onpe0() )
        cout << " randomize_v: temperature must be positive" << endl;
      return 1;
    }

    s->atoms.randomize_velocities(temp);
    return 0;
  }
};
#endif
