////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2016 The Regents of the University of California
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
// ResponseCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RESPONSECMD_H
#define RESPONSECMD_H

#include <iostream>
#include "UserInterface.h"

class Sample;
class ResponseCmd : public Cmd
{
  private:

  public:

  Sample *s;

  ResponseCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "response"; }
  const char *help_msg(void) const
  {
    return
    "\n response\n\n"
    " syntax: response amplitude nitscf [nite]\n\n"
    " syntax: response -vext vext_file nitscf [nite]\n\n"
    "   The response command computes the first order change in dipole caused\n"
    "   by an external electric field defined by the amplitude argument.\n"
    "   If the -vext option is used, the response command computes the\n"
    "   response to the external potential defined in the file vext_file.\n\n";
  }

  int action(int argc, char **argv);

};
#endif
