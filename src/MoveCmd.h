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
// MoveCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MOVECMD_H
#define MOVECMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class MoveCmd : public Cmd
{
  public:

  Sample *s;

  MoveCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "move"; }
  const char *help_msg(void) const
  {
    return
    "\n move\n\n"
    " syntax: move atom_name {to|by} x y z \n\n"
    "   The move command displaces an atom to a new position.\n"
    "   The new position is defined by absolute coordinates (to) or\n"
    "   by a relative displacement (by).\n"
    "   When using 'to', if one or more of the arguments is '*',\n"
    "   the corresponding component of the position is unchanged.\n\n";
  }

  int action(int argc, char **argv)
  {
    // move must have 6 arguments including the command name
    if ( argc != 6 )
    {
      if ( ui->onpe0() )
        cout << " use: move atom_name {to|by} x y z " << endl;
      return 1;
    }

    const string atom_name = argv[1];
    const string mode = argv[2];
    const string xs = argv[3];
    const string ys = argv[4];
    const string zs = argv[5];

    Atom* pa = s->atoms.findAtom(atom_name);
    if ( !pa )
    {
      if ( ui->onpe0() )
        cout << " MoveCmd: could not find atom " << atom_name << endl;
      return 1;
    }

    D3vector pos = pa->position();

    double x = pos.x;
    double y = pos.y;
    double z = pos.z;

    if ( mode == "to" )
    {
      if ( xs != "*" )
        x = atof(argv[3]);
      if ( ys != "*" )
        y = atof(argv[4]);
      if ( zs != "*" )
        z = atof(argv[5]);
      pos = D3vector(x,y,z);
    }
    else if ( mode == "by" )
    {
      x += atof(argv[3]);
      y += atof(argv[4]);
      z += atof(argv[5]);
      pos = D3vector(x,y,z);
    }
    else
    {
      if ( ui->onpe0() )
        cout << " MoveCmd: unknown mode" << endl;
      return 1;
    }

    pa->set_position(pos);
    if ( ui->onpe0() )
      cout << " MoveCmd: atom " << atom_name << " moved to "
           << pos << endl;

    return 0;
  }
};
#endif
