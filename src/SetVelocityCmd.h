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
// SetVelocityCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SETVELOCITYCMD_H
#define SETVELOCITYCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SetVelocityCmd : public Cmd
{
  public:

  Sample *s;

  SetVelocityCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "set_velocity"; }
  const char *help_msg(void) const
  {
    return
    "\n set_velocity\n\n"
    " syntax: set_velocity atom_name {vx|*|-}  {vy|*|-}  {vz|*|-}\n\n"
    "   The set_velocity command sets the velocity of an atom to vx,vy,vz.\n"
    "   If an argument is '*', the corresponding component of\n"
    "   the velocity is unchanged.\n"
    "   If an argument is '-', the corresponding component of\n"
    "   the velocity changes sign.\n\n";
  }

  int action(int argc, char **argv)
  {
    // set_velocity must have 5 arguments including the command name
    if ( argc != 5 )
    {
      if ( ui->onpe0() )
        cout << " use: set_velocity atom_name vx vy vz " << endl;
      return 1;
    }

    const string atom_name = argv[1];
    const string vxs = argv[2];
    const string vys = argv[3];
    const string vzs = argv[4];

    Atom* pa = s->atoms.findAtom(atom_name);
    if ( !pa )
    {
      if ( ui->onpe0() )
        cout << " SetVelocityCmd: could not find atom " << atom_name << endl;
      return 1;
    }

    D3vector vel = pa->velocity();

    double vx = vel.x;
    double vy = vel.y;
    double vz = vel.z;

    // change component only if argument is not "*"
    if ( vxs != "*" )
      vx = atof(argv[2]);
    if ( vys != "*" )
      vy = atof(argv[3]);
    if ( vzs != "*" )
      vz = atof(argv[4]);

    // change sign of component if argument is "-"
    if ( vxs == "-" )
      vx = -vel.x;
    if ( vys == "-" )
      vy = -vel.y;
    if ( vzs == "-" )
      vz = -vel.z;

    vel = D3vector(vx,vy,vz);

    pa->set_velocity(vel);
    if ( ui->onpe0() )
      cout << " SetVelocityCmd: atom " << atom_name << " set to "
           << vel << endl;

    return 0;
  }
};
#endif
