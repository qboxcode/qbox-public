////////////////////////////////////////////////////////////////////////////////
//
// SpeciesCmd.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesCmd.h,v 1.2 2007-10-19 16:24:05 fgygi Exp $

#ifndef SPECIESCMD_H
#define SPECIESCMD_H

#include <iostream>
#include <cstdlib>
#include <string>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class SpeciesCmd : public Cmd
{
  public:

  Sample *s;

  SpeciesCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "species"; }

  char *help_msg(void) const
  {
    return
    "\n species\n\n"
    " syntax: species name uri\n\n"
    "   The species command defines a species name.\n\n";
  }

  int action(int argc, char **argv);
};
#endif
