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

  void responseEfield(double amplitude, bool rpa, bool ipa,
    int nitscf, int nite);
  void responseVext(bool rpa, bool ipa, int nitscf, int nite, std::string fmt);

  public:

  Sample *s;

  ResponseCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "response"; }
  const char *help_msg(void) const
  {
    return
    "\n response\n\n"
    " syntax: response amplitude [-RPA|-IPA] nitscf [nite]\n"
    "         response -vext vext_file [-cube] [-RPA|-IPA] [-amplitude a] \n"
    "                  nitscf [nite]\n\n"
    "   The response command computes the polarizability tensor by\n"
    "   finite differences using external electric fields in the x,y,z\n"
    "   directions with magnitude defined by the amplitude argument.\n"
    "   If the -vext option is used, the response command computes the\n"
    "   response to the external potential defined in the file vext_file.\n"
    "   and writes the density response on a file vext_file.response.\n"
    "   Options:\n"
    "   -cube Use cube file format instead of XML for vext and response.\n"
    "   -RPA  Compute response within the Random Phase Approximation,\n"
    "         Vxc is frozen.\n"
    "   -IPA  Compute response within the Independent Particle Approximation,\n"
    "         VHartree and Vxc are frozen.\n"
    "   -amplitude a\n"
    "         Multiply the external potential by a before any calculation, \n"
    "         then divide the charge density response by a before output.\n\n";
  }

  int action(int argc, char **argv);

};
#endif
