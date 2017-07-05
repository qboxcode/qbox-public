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

#include<iostream>
#include "UserInterface.h"

class Sample;
class ResponseCmd : public Cmd
{
  private:

  void responseVext(bool rpa, int nitscf, int nite,
                    string io, int nx, int ny, int nz);
  void responseEfield(double amplitude, int nitscf, int nite);

  public:

  Sample *s;

  ResponseCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "response"; }
  const char *help_msg(void) const
  {
    return
    "\n response\n\n"
    " syntax: response amplitude nitscf [nite]\n"
    "         response -vext vext_file [-RPA] [-amplitude a] \n"
    "                  [-io iomode -nx nx -ny ny -nz nz] nitscf [nite]\n\n"
    "   The response command computes the polarizability tensor by\n"
    "   finite differences using external electric fields in the x,y,z\n"
    "   directions with magnitude defined by the amplitude argument.\n"
    "   If the -vext option is used, the response command computes the\n"
    "   response to the external potential defined in the file vext_file.\n"
    "   Acceptable control flags are:\n"
    "     1. -RPA. Compute response within Random Phase Approximation, Vxc is frozen.\n"
    "     2. -amplitude a. Scale the Vext by a before any calculations, \n"
    "                      then scale the charge density response by 1/a before output.\n"
    "     3. -io iomode. How the vext/response file shall be read/write. Possible choice:\n"
    "        cube: Gaussian cube format.\n"
    "        base64_serial or base64_parallel: base64-encoded binary grid function.\n"
    "                                          grid size need to be specified by nx, ny, nz.\n"
    "        use cube for best compatibility and base64_parallel for best performance\n\n";
  }

  int action(int argc, char **argv);

};
#endif
