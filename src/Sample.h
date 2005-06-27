////////////////////////////////////////////////////////////////////////////////
//
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Sample.h,v 1.7 2005-06-27 22:19:31 fgygi Exp $

#ifndef SAMPLE_H
#define SAMPLE_H

#include "AtomSet.h"
#include "ConstraintSet.h"
#include "Wavefunction.h"
#include "Control.h"

class Context;

class Sample
{
  private:
  
  public:
  
  const Context& ctxt_;

  AtomSet atoms;
  ConstraintSet constraints;
  Wavefunction wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;

  Sample(const Context& ctxt) : ctxt_(ctxt), atoms(ctxt), constraints(ctxt),
    wf(ctxt), wfv(0)
  { ctrl.sigmas = 0.5; ctrl.facs = 2.0; }
};
#endif
