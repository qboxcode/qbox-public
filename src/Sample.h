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
// Sample.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SAMPLE_H
#define SAMPLE_H

#include "AtomSet.h"
#include "ConstraintSet.h"
#include "ExtForceSet.h"
#include "Wavefunction.h"
#include "Control.h"

class Context;
class UserInterface;

class Sample
{
  private:

  public:

  const Context& ctxt_;

  AtomSet atoms;
  ConstraintSet constraints;
  ExtForceSet extforces;
  Wavefunction wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;
  UserInterface *ui;

  Sample(const Context& ctxt, UserInterface *ui_) : ctxt_(ctxt), ui(ui_),
    atoms(ctxt), constraints(ctxt),
    extforces(ctxt), wf(ctxt), wfv(0) {}
  ~Sample(void) { delete wfv; }
  void reset(void)
  {
    atoms.reset();
    constraints.reset();
    extforces.reset();
    wf.reset();
    delete wfv;
  }
};
#endif
