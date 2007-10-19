////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.C,v 1.6 2007-10-19 16:24:04 fgygi Exp $

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include "BOSampleStepper.h"
#include "CPSampleStepper.h"

#include<ctime>
#include<cassert>

int RunCmd::action(int argc, char **argv)
{

  if ( argc < 2 || argc > 4)
  {
    if ( ui->onpe0() )
      cout << " use: run niter" << endl;
      cout << "      run niter nitscf" << endl;
      cout << "      run niter nitscf nite" << endl;
    return 1;
  }

  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " <!-- RunCmd: no states, cannot run -->" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " <!-- RunCmd: ecut = 0.0, cannot run -->" << endl;
    return 1;
  }

  SampleStepper* stepper;

  int niter = atoi(argv[1]);
  int nite = 1;
  int nitscf = 1;
  if ( argc == 3 )
  {
    // run niter nitscf
    nitscf = atoi(argv[2]);
  }
  else if ( argc == 4 )
  {
    // run niter nitscf nite
    nitscf = atoi(argv[2]);
    nite = atoi(argv[3]);
  }
  if ( s->ctrl.wf_dyn == "MD" )
    stepper = new CPSampleStepper(*s);
  else
    stepper = new BOSampleStepper(*s,nitscf,nite);

  assert(stepper!=0);

  s->wf.info(cout,"wavefunction");
  stepper->step(niter);

  delete stepper;

  return 0;
}
