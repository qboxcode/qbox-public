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

#include "MPIdata.h"
#include "Context.h"
#include "AtomSet.h"
#include "ConstraintSet.h"
#include "ExtForceSet.h"
#include "Wavefunction.h"
#include "Control.h"

class UserInterface;
class ExternalPotential;

class Sample
{
  private:

  public:

  Context sd_ctxt;
  AtomSet atoms;
  ConstraintSet constraints;
  ExtForceSet extforces;
  ExternalPotential* vext;
  Wavefunction wf;
  Wavefunction* wfv; // wavefunction velocity
  Control ctrl;
  UserInterface *ui;

  Sample(UserInterface *ui_ = 0) : ui(ui_),
    sd_ctxt(MPIdata::sd_comm(),MPIdata::ngb(),MPIdata::nstb()),
    wf(sd_ctxt), wfv(0), vext(0) {}
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
