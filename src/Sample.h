////////////////////////////////////////////////////////////////////////////////
//
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Sample.h,v 1.6 2004-02-04 19:55:16 fgygi Exp $

#ifndef SAMPLE_H
#define SAMPLE_H

#include "AtomSet.h"
#include "Wavefunction.h"
#include "Control.h"

class Context;

class Sample
{
  private:
  
  public:
  
  const Context& ctxt_;

  AtomSet atoms;
  Wavefunction wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;

  Sample(const Context& ctxt) : ctxt_(ctxt), atoms(ctxt), wf(ctxt), wfv(0)
  { ctrl.sigmas = 0.5; ctrl.facs = 2.0; }
};
#endif
