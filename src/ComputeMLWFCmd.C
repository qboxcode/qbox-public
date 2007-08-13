////////////////////////////////////////////////////////////////////////////////
//
// ComputeMLWFCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ComputeMLWFCmd.C,v 1.1 2007-08-13 21:26:27 fgygi Exp $

#include "ComputeMLWFCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
using namespace std;

int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;
  SlaterDet& sd = *(wf.sd(0,0));
  
  mlwft = new MLWFTransform(sd);
  
  mlwft->compute_transform();
  mlwft->apply_transform(sd);
}
