////////////////////////////////////////////////////////////////////////////////
//
// QuitCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: QuitCmd.h,v 1.1 2003-01-10 00:35:36 fgygi Exp $

#ifndef QUITCMD_H
#define QUITCMD_H

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>

#include "UserInterface.h"
#include "Sample.h"

class QuitCmd : public Cmd
{
  public:

  Sample* s;

  QuitCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "quit"; }

  char *help_msg(void) const
  {
    return 
    "\n quit\n\n"
    " syntax: quit\n\n"
    "   The quit command exits without saving any data.\n\n";
  }

  int action(int argc, char **argv)
  {      
    ui->terminate();
    return 0;
  }
};
#endif
