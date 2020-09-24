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
// KpointCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef KPOINTCMD_H
#define KPOINTCMD_H

#include "MPIdata.h"
#include "UserInterface.h"
#include "D3vector.h"
#include "Sample.h"

class KpointCmd : public Cmd
{
  public:

  Sample *s;

  KpointCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "kpoint"; }
  const char *help_msg(void) const
  {
    return
    "\n kpoint\n\n"
    " syntax:\n\n"
    "   kpoint add kx ky kz weight\n"
    "   kpoint delete kx ky kz \n"
    "   kpoint move kx ky kz new_kx new_ky new_kz\n"
    "   kpoint list\n\n";
  }

  int action(int argc, char **argv)
  {
    const bool onpe0 = MPIdata::onpe0();
    if ( argc < 2 )
    {
      if ( onpe0 )
        cout << help_msg();
      return 1;
    }
    string subcmd(argv[1]);
    if ( subcmd == "add" )
    {
      if ( argc != 6 )
      {
        if ( onpe0 )
          cout << help_msg();
        return 1;
      }
      double kx = atof(argv[2]);
      double ky = atof(argv[3]);
      double kz = atof(argv[4]);
      double w = atof(argv[5]);
      s->wf.add_kpoint(D3vector(kx,ky,kz),w);
    }
    else if ( subcmd == "delete" )
    {
      if ( argc != 5 )
      {
        if ( onpe0 )
          cout << help_msg();
        return 1;
      }
      double kx = atof(argv[2]);
      double ky = atof(argv[3]);
      double kz = atof(argv[4]);
      s->wf.del_kpoint(D3vector(kx,ky,kz));
    }
    else if ( subcmd == "move" )
    {
      if ( argc != 8 )
      {
        if ( onpe0 )
          cout << help_msg();
        return 1;
      }
      double kx = atof(argv[2]);
      double ky = atof(argv[3]);
      double kz = atof(argv[4]);
      double newkx = atof(argv[5]);
      double newky = atof(argv[6]);
      double newkz = atof(argv[7]);
      s->wf.move_kpoint(D3vector(kx,ky,kz),D3vector(newkx,newky,newkz));
    }
    else if ( subcmd == "list" )
    {
      if ( argc != 2 )
      {
        if ( onpe0 )
          cout << help_msg();
        return 1;
      }
      if ( onpe0 )
      {
        cout << " kpoint list: reciprocal lattice coordinates" << endl;
        for ( int ikp = 0; ikp < s->wf.nkp(); ikp++ )
        {
          cout << " "
               << s->wf.kpoint(ikp) << "    " << s->wf.weight(ikp) << endl;
        }
        cout << " kpoint list: cartesian coordinates" << endl;
        const UnitCell& u = s->wf.cell();
        double wsum = 0.0;
        for ( int ikp = 0; ikp < s->wf.nkp(); ikp++ )
        {
          D3vector kp = s->wf.kpoint(ikp);
          double w = s->wf.weight(ikp);
          cout << " "
               << kp * u.b(0) << " " << kp * u.b(1) << " " << kp * u.b(2)
               << "    " << w << endl;
          wsum += w;
        }
        cout << " total weight: " << wsum << endl;
      }
    }
    else
    {
      if ( onpe0 )
        cout << help_msg();
    }

    return 0;
  }
};
#endif
