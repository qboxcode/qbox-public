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
// RseedCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RSEEDCMD_H
#define RSEEDCMD_H

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>

#include "UserInterface.h"
#include "Sample.h"

class RseedCmd : public Cmd
{
  public:

  Sample* s;

  RseedCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "rseed"; }
  const char *help_msg(void) const
  {
    return
    "\n rseed\n\n"
    " syntax: rseed [seed_value]\n\n"
    "   The rseed command initializes the random number generator.\n"
    "   If no argument is given, the time() function is used as a seed value"
    "\n\n";
  }

  int action(int argc, char **argv)
  {
    int seed = (int) time(0);
    if ( argc == 2 )
      seed = atoi(argv[1]);
    else
    {
      if ( ui->onpe0() )
      {
        s->ctxt_.ibcast_send(1,1,&seed,1);
        cout << "<seed> " << seed << " </seed>" << endl;
      }
      else
      {
        s->ctxt_.ibcast_recv(1,1,&seed,1,0,0);
      }
    }
    srand48((long int)seed);
    return 0;
  }
};
#endif
