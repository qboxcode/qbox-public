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
// StatusCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STATUSCMD_H
#define STATUSCMD_H

#include <iostream>
#include "UserInterface.h"
#include "Sample.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"

class StatusCmd : public Cmd
{
  private:

  int niter, nfi;

  public:

  Sample *s;

  StatusCmd(Sample *sample) : s(sample) {};

  const char *name(void) const { return "status"; }
  const char *help_msg(void) const
  {
    return
    "\n run\n\n"
    " syntax: status \n\n"
    "   The status command print information about the current\n"
    "   status of the simulation\n\n";
  }

  int action(int argc, char **argv)
  {
    // compute the size of the potential grid
    int np0v = 0;
    int np1v = 0;
    int np2v = 0;
    if ( s->wf.ecut() > 0 && s->wf.cell().volume() > 0 )
    {
      ChargeDensity cd(s->wf);
      np0v = cd.vft()->np0();
      np1v = cd.vft()->np1();
      np2v = cd.vft()->np2();
    }
    if ( ui->onpe0() )
    {
      cout << "<np0v> " << np0v << " </np0v>  "
           << "<np1v> " << np1v << " </np1v>  "
           << "<np2v> " << np2v << " </np2v>" << endl;
    }
    s->wf.info(cout,"wf");
    if ( s->wfv != 0 )
      s->wfv->info(cout,"wfv");
    if ( ui->onpe0() )
    {
      setprecision(8);
      cout << "<rcm> " << setprecision(8) << s->atoms.rcm()
           << " </rcm>" << endl;
      cout << "<vcm> " << setprecision(8) << s->atoms.vcm()
           << " </vcm>" << endl;
    }
    return 0;
  }
};
#endif
