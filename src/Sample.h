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
#include "AtomSet.h"
#if 0
#include "ConstraintSet.h"
#include "ExtForceSet.h"
#include "Wavefunction.h"
#endif
#include "Control.h"

class UserInterface;
class ExternalPotential;

class Sample
{
  private:

  public:

  AtomSet atoms;
  //ConstraintSet constraints;
  //ExtForceSet extforces;
  ExternalPotential* vext;
  //Wavefunction wf;
  //Wavefunction* wfv; // wavefunction velocity
  Control ctrl;
  UserInterface *ui;

  //Sample(UserInterface *ui_ = 0) : ui(ui_), atoms(ctxt), constraints(ctxt),
  //  extforces(ctxt), vext(0), wf(ctxt), wfv(0) {}
  Sample(UserInterface *ui_ = 0) : ui(ui_) {}
  //~Sample(void) { delete wfv; }
  void reset(void)
  {
    atoms.reset();
    //constraints.reset();
    //extforces.reset();
    //wf.reset();
    //delete wfv;
  }
};
#endif
