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
// RunCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include "BOSampleStepper.h"
#include "CPSampleStepper.h"

#include<ctime>
#include<cassert>
#include<string>
#include<stdexcept>

int RunCmd::action(int argc, char **argv)
{
  string usage(" use: run [-atomic_density] niter"
               "      run [-atomic_density] niter nitscf"
               "      run [-atomic_density] niter nitscf nite");
  if ( argc < 2 || argc > 5)
    throw invalid_argument(usage);
  if ( s->wf.nst() == 0 )
    throw runtime_error("RunCmd: no states, cannot run");
  if ( s->wf.ecut() == 0.0 )
    throw runtime_error("RunCmd: ecut = 0, cannot run");
  if ( s->wf.cell().volume() == 0.0 )
    throw runtime_error("RunCmd: volume = 0, cannot run");

  SampleStepper* stepper;

  int iarg = 1;
  bool atomic_density = false;
  if ( !strcmp(argv[iarg],"-atomic_density") )
  {
    atomic_density = true;
    iarg++;
    argc--;
  }

  int niter = atoi(argv[iarg]);
  int nite = 0;
  int nitscf = 1;
  if ( argc == 3 )
  {
    // run niter nitscf
    nitscf = atoi(argv[iarg+1]);
  }
  else if ( argc == 4 )
  {
    // run niter nitscf nite
    nitscf = atoi(argv[iarg+1]);
    nite = atoi(argv[iarg+2]);
  }

  s->extforces.setup(s->atoms);

  if ( s->ctrl.wf_dyn == "MD" )
    stepper = new CPSampleStepper(*s);
  else
    stepper = new BOSampleStepper(*s,nitscf,nite);

  assert(stepper!=0);
  stepper->set_iter_cmd(s->ctrl.iter_cmd);
  stepper->set_iter_cmd_period(s->ctrl.iter_cmd_period);

  if ( atomic_density )
    stepper->initialize_density();

  s->wf.info(cout,"wavefunction");
  stepper->step(niter);

  // Delete wave function velocity if not using atoms_dyn = MD
  if ( s->ctrl.atoms_dyn != "MD" )
  {
    if ( s->wfv != 0 )
      delete s->wfv;
    s->wfv = 0;
  }

  delete stepper;

  return 0;
}
