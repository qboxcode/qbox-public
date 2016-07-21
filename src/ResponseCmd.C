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
// ResponseCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include "ResponseCmd.h"
#include<iostream>
using namespace std;
#include "BOSampleStepper.h"

#include<ctime>
#include<cassert>

int ResponseCmd::action(int argc, char **argv)
{
  if ( argc < 3 || argc > 4)
  {
    if ( ui->onpe0() )
      cout << " use: response vext_file nitscf [nite]" << endl;
    return 1;
  }

  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " ResponseCmd: no states, cannot run" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " ResponseCmd: ecut = 0.0, cannot run" << endl;
    return 1;
  }

  string vext_filename(argv[1]);

  int niter = atoi(argv[2]);
  int nite = 0;
  if ( argc == 4 )
  {
    // nitscf nite
    nite = atoi(argv[3]);
  }

  SampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);
  assert(stepper!=0);

  // update external potential from file
  if ( !s.vext )
  {
    if ( ui->onpe0() )
      cout << " ResponseCmd: vext file not defined" << endl;
    return 1;
  }

  ifstream vextfile(s.vext.c_str());

  stepper->step(0);
  // compute density

  // change sign of potential

  stepper->step(0);
  // compute density

  // compute density difference
  // output density difference

  delete stepper;

  return 0;
}
