////////////////////////////////////////////////////////////////////////////////
//
// SaveCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveCmd.h,v 1.3 2008-01-26 01:34:11 fgygi Exp $

#ifndef SAVECMD_H
#define SAVECMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SaveCmd : public Cmd
{
  public:

  Sample *s;

  SaveCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "save"; }

  char *help_msg(void) const
  {
    return
    "\n save\n\n"
    " syntax: save [-serial] [-text] filename \n\n"
    "   The save command saves the sample to the file filename.\n\n"
    "   When using the -serial option, I/O is performed from the \n"
    "   head node only. \n\n";
  }

  int action(int argc, char **argv);
};
#endif
