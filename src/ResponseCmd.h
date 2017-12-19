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

  void responseEfield(double amplitude, int nitscf, int nite);
  void responseVext(bool rpa, bool ipa, int nitscf, int nite, string io);

  public:

  Sample *s;

  ResponseCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "response"; }
  const char *help_msg(void) const
  {
    return
    "\n response\n\n"
    " syntax: response amplitude nitscf [nite]\n"
    "         response -vext vext_file [-RPA|-IPA] [-amplitude a] \n"
    "                  [-io iomode -nx nx -ny ny -nz nz] nitscf [nite]\n\n"
    "   The response command computes the polarizability tensor by\n"
    "   finite differences using external electric fields in the x,y,z\n"
    "   directions with magnitude defined by the amplitude argument.\n"
    "   If the -vext option is used, the response command computes the\n"
    "   response to the external potential defined in the file vext_file.\n"
    "   Control flags:\n"
    "   -RPA  Compute response within the Random Phase Approximation,\n"
    "         Vxc is frozen.\n"
    "   -IPA  Compute response within the Independent Particle Approximation,\n"
    "         VHartree and Vxc are frozen.\n"
    "   -amplitude a\n"
    "         Scale the external potential by a before any calculation, \n"
    "         then scale the charge density response by 1/a before output.\n"
    "   -io iomode\n"
    "         Valid choices of iomode: cube, base64_serial, base64_parallel\n"
    "   -nx nx, -ny ny, -nz nz\n"
    "         grid size (for base64_serial and base64_parallel only)\n\n";
  }

  int action(int argc, char **argv);

};
#endif
