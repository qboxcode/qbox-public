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
// CPSampleStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "CPSampleStepper.h"
#include "SlaterDet.h"
#include "MDWavefunctionStepper.h"
#include "MDIonicStepper.h"
#include "SDCellStepper.h"
#include "CGCellStepper.h"
#include "Basis.h"
#include "Species.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CPSampleStepper::CPSampleStepper(Sample& s) :
  SampleStepper(s), cd_(s.wf), ef_(s,cd_), dwf(s.wf), wfv(s.wfv)
{
  const double emass = s.ctrl.emass;
  const double dt = s.ctrl.dt;
  double dt2bye = (emass == 0.0) ? 0.5 / s.wf.ecut() : dt*dt/emass;

  // divide dt2bye by facs coefficient if stress == ON
  const double facs = 2.0;
  if ( s.ctrl.stress == "ON" )
  {
    dt2bye /= facs;
  }
  if ( s.wfv == 0 )
  {
    s.wfv = new Wavefunction(s.wf);
    s.wfv->clear();
  }
  mdwf_stepper = new MDWavefunctionStepper(s.wf,s.wfv,dt,dt2bye,tmap);
  assert(mdwf_stepper!=0);
  mdionic_stepper = 0;
  if ( s.ctrl.atoms_dyn != "LOCKED" )
    mdionic_stepper = new MDIonicStepper(s);
}

////////////////////////////////////////////////////////////////////////////////
CPSampleStepper::~CPSampleStepper(void)
{
  delete mdwf_stepper;
  if ( mdionic_stepper != 0 ) delete mdionic_stepper;

  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = i->second.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( MPIdata::onpe0() && (tmax > 0.0) )
    {
      string s = "name=\"" + i->first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void CPSampleStepper::step(int niter)
{
  const bool onpe0 = MPIdata::onpe0();

  // check that there are no fractionally occupied states
  // next line: (3-nspin) = 2 if nspin==1 and 1 if nspin==2
  if ( s_.wf.nel() != (( 3 - s_.wf.nspin() ) * s_.wf.nst()) )
  {
    if ( onpe0 )
    {
      cout << " CPSampleStepper::step:"
              " fractionally occupied or empty states: cannot run" << endl;
    }
    return;
  }

  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;

  const double dt = s_.ctrl.dt;
  double ekin_ion=0.0,ekin_e;

  const string wf_dyn = s_.ctrl.wf_dyn;
  assert(wf_dyn=="MD");
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;

  const bool compute_hpsi = true;
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);
  else if ( cell_dyn == "CG" )
    cell_stepper = new CGCellStepper(s_);

  if ( s_.wfv == 0 )
  {
    s_.wfv = new Wavefunction(wf);
    s_.wfv->clear();
  }

  if ( mdionic_stepper )
    mdionic_stepper->setup_constraints();

  Timer tm_iter;

  tmap["charge"].start();
  cd_.update_density();
  tmap["charge"].stop();

  ef_.update_vhxc(compute_stress);
  double energy =
    ef_.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma_eks);

  mdwf_stepper->compute_wfm(dwf);

  for ( int iter = 0; iter < niter; iter++ )
  {
    tm_iter.reset();
    tm_iter.start();
    if ( onpe0 )
      cout << "<iteration count=\"" << iter+1 << "\">\n";

    mdwf_stepper->update(dwf);

    ekin_e = mdwf_stepper->ekin();

    if ( onpe0 )
    {
      cout << ef_;
      if ( ef_.el_enth() )
        cout << *ef_.el_enth();
    }

    if ( compute_forces )
    {
      if ( iter > 0 )
      {
        mdionic_stepper->compute_v(energy,fion);
      }
      mdionic_stepper->compute_r(energy,fion);
      ekin_ion = mdionic_stepper->ekin();

      if ( onpe0 )
      {
        cout << "<atomset>" << endl;
        cout << atoms.cell();
        for ( int is = 0; is < atoms.atom_list.size(); is++ )
        {
          int i = 0;
          for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
          {
            Atom* pa = atoms.atom_list[is][ia];
            cout << "  <atom name=\"" << pa->name() << "\""
                 << " species=\"" << pa->species()
                 << "\">\n"
                 << "    <position> "
                 << mdionic_stepper->r0(is,i) << " "
                 << mdionic_stepper->r0(is,i+1) << " "
                 << mdionic_stepper->r0(is,i+2) << " </position>\n"
                 << "    <velocity> "
                 << mdionic_stepper->v0(is,i) << " "
                 << mdionic_stepper->v0(is,i+1) << " "
                 << mdionic_stepper->v0(is,i+2) << " </velocity>\n"
                 << "    <force> "
                 << fion[is][i] << " "
                 << fion[is][i+1] << " "
                 << fion[is][i+2]
                 << " </force>\n  </atom>" << endl;

            i += 3;
          }
        }
        cout << "</atomset>" << endl;
      }

      if ( s_.constraints.size() > 0 )
      {
        s_.constraints.compute_forces(mdionic_stepper->r0(), fion);
        if ( onpe0 )
        {
          s_.constraints.list_constraints(cout);
        }
      }

    }

    if ( onpe0 )
    {
      cout << "  <ekin_e> " << ekin_e << " </ekin_e>\n";

      if ( compute_forces )
      {
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
        cout << "  <temp_ion> " << mdionic_stepper->temp() << " </temp_ion>\n";
        cout << "  <eta_ion> " << mdionic_stepper->eta() << " </eta_ion>\n";
      }
      double econst = energy + ekin_ion + ekin_e;
      if ( mdionic_stepper )
        econst += mdionic_stepper->ekin_stepper();
      cout << "  <econst> " << econst << " </econst>\n";
      cout << "  <ekin_ec> " << econst + ekin_e << " </ekin_ec>\n";
    }

    if ( compute_stress )
    {
      compute_sigma();
      print_stress();

      if ( cell_dyn != "LOCKED" )
      {
        cell_stepper->compute_new_cell(energy,sigma,fion);

        // Update cell
        cell_stepper->update_cell();

        ef_.cell_moved();
        ef_.atoms_moved();
      }
    }

    if ( compute_forces )
    {
      ef_.atoms_moved();
    }

    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    ef_.update_vhxc(compute_stress);
    energy =
      ef_.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma_eks);
    ef_.enthalpy();

    tm_iter.stop();
    // print iteration time
    double time = tm_iter.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( onpe0 )
    {
      string s = "name=\"iteration\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }

    if ( onpe0 )
      cout << "</iteration>" << endl;

    if ( compute_forces )
      s_.constraints.update_constraints(dt);
  } // iter

  // dwf and fion now contain the forces on wavefunctions and ions at the
  // endpoint

  mdwf_stepper->compute_wfv(dwf); // replace wfm by wfv

  if ( compute_forces )
  {
    // Note: next line function call updates velocities in the AtomSet
    mdionic_stepper->compute_v(energy,fion);
  }

  if ( cell_stepper != 0 ) delete cell_stepper;
}
