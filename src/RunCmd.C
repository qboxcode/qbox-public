////////////////////////////////////////////////////////////////////////////////
//
// RunCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RunCmd.C,v 1.2 2003-02-04 19:21:30 fgygi Exp $

#include "RunCmd.h"
#include<iostream>
using namespace std;
#include "EnergyFunctional.h"
#include "SampleStepper.h"

#include<ctime>
#include<cassert>

int RunCmd::action(int argc, char **argv)
{

  if ( argc != 2 )
  {
    if ( ui->onpe0() )
      cout << " use: run nsteps" << endl;
    return 1;
  }
  
  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " <qb:error> RunCmd: no states, cannot run </qb:error>" << endl;
    return 1;
  }
  
  EnergyFunctional ef(*s);
  SampleStepper stepper(*s);
  
  int niter = atoi(argv[1]);
  stepper.step(ef,niter);
  
  return 0;
}
