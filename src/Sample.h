////////////////////////////////////////////////////////////////////////////////
//
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Sample.h,v 1.5 2003-06-11 22:10:11 fgygi Exp $

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

  Sample(const Context& ctxt) : ctxt_(ctxt), atoms(ctxt), wf(ctxt), wfv(0) {}
};
#endif
