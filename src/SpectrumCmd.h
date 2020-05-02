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
// SpectrumCmd.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SPECTRUMCMD_H
#define SPECTRUMCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
class Sample;

class SpectrumCmd : public Cmd
{
  private:

  public:
  Sample *s;

  SpectrumCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "spectrum"; }
  const char *help_msg(void) const
  {
    return
    "\n spectrum\n\n"
    " syntax: spectrum filename\n"
    "         spectrum width filename\n"
    "         spectrum emin emax width filename\n\n"
    "   The spectrum command computes the dipole transition strengths\n"
    " between occupied and empty orbitals. It computes Kohn-Sham eigenvalues\n"
    " and eigenfunctions of the current wave function using the the current\n"
    " value of the xc variable. The corresponding absorption\n"
    " spectrum is written on an output file after convolution with a\n"
    " gaussian function.\n"
    "   emin, emax: energy range (optional)\n"
    "   width     : width of the gaussian (optional) (default 0.05 eV)\n"
    "   filename  : output file name\n"
    " If emin and emax are not given, the energy range includes\n"
    " all possible transitions between occupied and empty orbitals.\n"
    " All energy parameters must be given in eV.\n\n";
  }

  int action(int argc, char **argv);

  SpectrumCmd();
};
#endif
