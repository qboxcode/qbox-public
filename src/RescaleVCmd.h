////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2010 The Regents of the University of California
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
// RescaleVCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RESCALEVCMD_H
#define RESCALEVCMD_H

#include "UserInterface.h"
#include "Sample.h"
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

class RescaleVCmd : public Cmd
{
  private:

  public:

  Sample *s;

  RescaleVCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "rescale_v"; }
  const char *help_msg(void) const
  {
    return
    "\n rescale_v\n\n"
    " syntax: rescale_v f\n\n"
    "   The rescale_v command multiplies the velocity of all atoms \n"
    "   by the factor f.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
        cout << " use: rescale_v fac" << endl;
      return 1;
    }
    const double fac = atof(argv[1]);
    s->atoms.rescale_velocities(fac);
    return 0;
  }
};
#endif
