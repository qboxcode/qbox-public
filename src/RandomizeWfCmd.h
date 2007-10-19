////////////////////////////////////////////////////////////////////////////////
//
// RandomizeWfCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RandomizeWfCmd.h,v 1.3 2007-10-19 16:24:04 fgygi Exp $

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

  char *name(void) const { return "randomize_wf"; }
  char *help_msg(void) const
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
