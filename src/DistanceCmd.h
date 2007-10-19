////////////////////////////////////////////////////////////////////////////////
//
// DistanceCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: DistanceCmd.h,v 1.3 2007-10-19 16:24:04 fgygi Exp $

#ifndef DISTANCECMD_H
#define DISTANCECMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include <cstdlib>

class DistanceCmd : public Cmd
{
  public:

  Sample *s;

  DistanceCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "distance"; }
  char *help_msg(void) const
  {
    return
    "\n distance\n\n"
    " syntax: distance name1 name2\n\n"
    "   The distance command prints the distance between two atoms.\n\n";
  }

  int action(int argc, char **argv)
  {
    if ( argc != 3 )
    {
      if ( ui->onpe0() )
      {
        cout << " use: distance name1 name2" << endl;
      }
      return 1;
    }

    string name1(argv[1]);
    string name2(argv[2]);
    Atom* a1 = s->atoms.findAtom(name1);
    Atom* a2 = s->atoms.findAtom(name2);
    if ( a1 == 0 || a2 == 0 )
    {
      // either a1 or a2 was not found
      if ( ui->onpe0() )
      {
        if ( a1 == 0 )
          cout << " <!-- DistanceCmd: atom " << name1 << " not defined -->"
               << endl;
        if ( a2 == 0 )
          cout << " <!-- DistanceCmd: atom " << name2 << " not defined -->"
               << endl;
      }
      return 1;
    }

    if ( ui->onpe0() )
    {
      const double d = length(a1->position()-a2->position());
      cout.setf(ios::fixed,ios::floatfield);
      cout << " <!-- distance " << name1 << "-" << name2 << ": "
           << setprecision(3)
           << d << " (a.u.) / " << 0.529177*d << " (Ang) -->" << endl;
    }
    return 0;
  }
};
#endif
