////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.C,v 1.4 2004-03-11 21:52:31 fgygi Exp $

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include "EnergyFunctional.h"
#include "BOSampleStepper.h"
#include "CPSampleStepper.h"

#include<ctime>
#include<cassert>

int RunCmd::action(int argc, char **argv)
{

  if ( argc < 2 || argc > 3)
  {
    if ( ui->onpe0() )
      cout << " use: run niter" << endl;
    return 1;
  }
  
  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " <qb:error> RunCmd: no states, cannot run </qb:error>" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " <qb:error> RunCmd: ecut = 0.0, cannot run </qb:error>" << endl;
    return 1;
  }
  
  EnergyFunctional ef(*s);
  SampleStepper* stepper;
  
  int niter = atoi(argv[1]);
  int nite = 1;
  if ( argc == 3 )
  {
    // run niter nite
    nite = atoi(argv[2]);
  }
  
  if ( s->ctrl.wf_dyn == "MD" )
    stepper = new CPSampleStepper(*s,ef);
  else
    stepper = new BOSampleStepper(*s,ef,nite);
  
  assert(stepper!=0);
  
  s->wf.info(cout,"wavefunction");
  stepper->step(niter);
  
  delete stepper;
  
  return 0;
}
