////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "SDIonicStepper.h"
#include "MDIonicStepper.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nite) : SampleStepper(s), 
  dwf(s.wf), wfv(s.wfv), nite_(nite)
{
  fion.resize(s_.atoms.nsp());
  for ( int is = 0; is < fion.size(); is++ )
    fion[is].resize(3*s_.atoms.na(is));
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(EnergyFunctional& e, int niter)
{
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  UnitCell dcell;
  atoms.get_positions(tau0);
  atoms.get_velocities(vel);
  
  const double dt = s_.ctrl.dt;
  double ekin_ion=0.0, temp_ion=0.0, eta;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  
  const bool compute_hpsi = ( wf_dyn != "LOCKED" );
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );

  const bool ttherm = ( s_.ctrl.thermostat == "ON" );
  const int ndofs = 3 * s_.atoms.size();
  const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
  const double th_temp = s_.ctrl.th_temp;
  const double th_time = s_.ctrl.th_time;
  
  Timer tm_iter;
  
  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
    wf_stepper = new SDWavefunctionStepper(s_,tmap);
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(s_,tmap);
  assert(wf_stepper!=0);
    
  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "MD" )
    ionic_stepper = new MDIonicStepper(s_);
  
  for ( int iter = 0; iter < niter; iter++ )
  {
    // ionic iteration
 
    tm_iter.start();
 
    if ( s_.ctxt_.onpe0() )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
 
    // wavefunction extrapolation
    if ( s_.wfv != 0 && atoms_dyn != "LOCKED" )
    {
      //wf_stepper->extrapolate(wf,wfv);
    }

    double energy;
    // do nite-1 electronic steps without forces
    if ( wf_stepper != 0 )
    {
      for ( int ite = 0; ite < nite_-1; ite++ )
      {
        energy = e.energy(true,dwf,false,fion,false,dcell);
 
        wf_stepper->update(dwf);
 
        if ( s_.ctxt_.onpe0() )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <ekin>   " << setprecision(8)
               << setw(15) << e.ekin() << " </ekin>\n"
               << "  <eps>    " << setw(15) << e.eps() << " </eps>\n"
               << "  <enl>    " << setw(15) << e.enl() << " </enl>\n"
               << "  <ecoul>  " << setw(15) << e.ecoul() << " </ecoul>\n"
               << "  <exc>    " << setw(15) << e.exc() << " </exc>\n"
               << "  <esr>    " << setw(15) << e.esr() << " </esr>\n"
               << "  <eself>  " << setw(15) << e.eself() << " </eself>\n"
               << "  <etotal> " << setw(15) << e.etotal() << " </etotal>\n"
               << flush;
        }
      }
    }
 
    // last electronic iteration
    energy =
      e.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,dcell);

    if ( s_.ctxt_.onpe0() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout.setf(ios::right,ios::adjustfield);
      cout << "  <ekin>   " << setprecision(8)
           << setw(15) << e.ekin() << " </ekin>\n"
           << "  <eps>    " << setw(15) << e.eps() << " </eps>\n"
           << "  <enl>    " << setw(15) << e.enl() << " </enl>\n"
           << "  <ecoul>  " << setw(15) << e.ecoul() << " </ecoul>\n"
           << "  <exc>    " << setw(15) << e.exc() << " </exc>\n"
           << "  <esr>    " << setw(15) << e.esr() << " </esr>\n"
           << "  <eself>  " << setw(15) << e.eself() << " </eself>\n"
           << "  <etotal> " << setw(15) << e.etotal() << " </etotal>\n"
           << flush;
    }
 
    // make ionic step
    
    if ( compute_forces )
    {
      if ( s_.wf.context().onpe0() )
      {
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
                 << ionic_stepper->tau0(is,i) << " "
                 << ionic_stepper->tau0(is,i+1) << " " 
                 << ionic_stepper->tau0(is,i+2) << " </position>\n"
                 << "    <force> " << fion[is][i] << " "
                 << fion[is][i+1] << " " << fion[is][i+2]
                 << " </force>\n  </atom>" << endl;
 
            i += 3;
          }
        }
      }
 
      if ( iter == 0 )
        ionic_stepper->preprocess(fion);
 
      ionic_stepper->update(fion);
      ekin_ion = ionic_stepper->ekin();
 
      if ( iter == niter-1 )
        ionic_stepper->postprocess(fion);
        
      e.atoms_moved();
    }

    if ( s_.ctxt_.onpe0() )
    {
      cout << "  <econst> " << energy+ekin_ion << " </econst>\n";
    }
 
    // print iteration time
    double time = tm_iter.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.onpe0() )
    {
      cout << "  <!-- timing "
           << setw(15) << "iteration"
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
      cout << "  </iteration>" << endl;
    }
  } // for iter
}
