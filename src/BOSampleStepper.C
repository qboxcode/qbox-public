////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.2 2003-11-27 01:20:21 fgygi Exp $

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
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
  const double dt_inv = 1.0 / dt;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  
  const bool compute_hpsi = ( wf_dyn != "LOCKED" );
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );

  Timer tm_iter;
  
  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
    wf_stepper = new SDWavefunctionStepper(s_,tmap);
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(s_,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(s_,tmap);
  assert(wf_stepper!=0);
    
  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "MD" )
    ionic_stepper = new MDIonicStepper(s_);
  
  // Allocate wavefunction velocity if not available
  if ( atoms_dyn != "LOCKED" )
  {
    if ( s_.wfv == 0 )
    {
      s_.wfv = new Wavefunction(wf);
      s_.wfv->clear();
    }
  }
      
  for ( int iter = 0; iter < niter; iter++ )
  {
    // ionic iteration
 
    tm_iter.start();
 
    if ( s_.ctxt_.onpe0() )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
 
    if ( compute_forces )
    {
      // compute energy and ionic forces using existing wavefunction
      double energy =
        e.energy(false,dwf,compute_forces,fion,compute_stress,dcell);

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
 
      if ( iter == 0 )
        ionic_stepper->preprocess(fion);
        
      // save a copy of atomic positions at time t0
      vector<vector<double> > tau0tmp = ionic_stepper->tau0();
      
      ionic_stepper->update(fion);
      
      // ekin_ion is the kinetic energy at time t0
      const double ekin_ion = ionic_stepper->ekin();
      const double temp_ion = ionic_stepper->temp();
        
      // positions, velocities and forces at time t0 are now known
      // print velocities and forces at time t0
      if ( s_.ctxt_.onpe0() )
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
                 << "    <velocity> " 
                 << ionic_stepper->vel(is,i) << " "
                 << ionic_stepper->vel(is,i+1) << " " 
                 << ionic_stepper->vel(is,i+2) << " </velocity>\n"
                 << "    <force> " 
                 << fion[is][i] << " "
                 << fion[is][i+1] << " " 
                 << fion[is][i+2]
                 << " </force>\n  </atom>" << endl;
 
            i += 3;
          }
        }
        cout << "  <econst> " << energy+ekin_ion << " </econst>\n";
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
        cout << "  <temp_ion> " << ionic_stepper->temp() << " </temp_ion>\n";
      }
      
      e.atoms_moved();
    }
    
    // wavefunction extrapolation
    if ( compute_forces )
    {
      // extrapolate wavefunctions
      // s_.wfv contains the wavefunction velocity
      for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
      {
        for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
        {
          if ( s_.wf.sd(ispin,ikp) != 0 )
          {
            if ( s_.wf.sdcontext(ispin,ikp)->active() )
            {
              if ( iter == 0 )
              {
                // extrapolation using the wavefunction velocity
                // v = c (save current c in v)
                // c = c + dt * v
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                for ( int i = 0; i < 2*mloc*nloc; i++ )
                {
                  const double x = c[i];
                  const double v = cv[i];
                  c[i] = x + dt * v;
                  cv[i] = x;
                }
                tmap["gram"].start();
                s_.wf.sd(ispin,ikp)->gram();
                tmap["gram"].stop();
              }
              else
              {
                // extrapolation using the previous wavefunction
                // cm is stored in wfv
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                for ( int i = 0; i < 2*mloc*nloc; i++ )
                {
                  const double x = c[i];
                  const double xm = cv[i];
                  c[i] = 2.0 * x - xm;
                  cv[i] = x;
                }
                tmap["riccati"].start();
                s_.wf.sd(ispin,ikp)->riccati(*s_.wfv->sd(ispin,ikp));
                tmap["riccati"].stop();
              }
            }
          }
        }
      }
    }

    // do nite electronic steps without forces
    if ( wf_stepper != 0 )
    {
      wf_stepper->preprocess();
      for ( int ite = 0; ite < nite_; ite++ )
      {
        // at the last nite iteration, compute ionic forces for the last
        double energy = e.energy(true,dwf,false,fion,false,dcell);
 
        wf_stepper->update(dwf);
 
        if ( s_.ctxt_.onpe0() )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <ekin>   " << setprecision(8)
               << setw(15) << e.ekin() << " </ekin>\n"
               << "  <eps_int>    " << setw(15) << e.eps() << " </eps_int>\n"
               << "  <enl_int>    " << setw(15) << e.enl() << " </enl_int>\n"
               << "  <ecoul_int>  " << setw(15) << e.ecoul() << " </ecoul_int>\n"
               << "  <exc_int>    " << setw(15) << e.exc() << " </exc_int>\n"
               << "  <esr_int>    " << setw(15) << e.esr() << " </esr_int>\n"
               << "  <eself_int>  " << setw(15) << e.eself() << " </eself_int>\n"
               << "  <etotal_int> " << setw(15) << e.etotal() << " </etotal_int>\n"
               << flush;
        }
      }
      wf_stepper->postprocess();
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
  
  if ( compute_forces )
  {
    // compute wavefunction velocity after last iteration
    // s_.wfv contains the previous wavefunction
    for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
      {
        if ( s_.wf.sd(ispin,ikp) != 0 )
        {
          if ( s_.wf.sdcontext(ispin,ikp)->active() )
          {
            double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
            double* cm = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
            const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
            const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
            for ( int i = 0; i < 2*mloc*nloc; i++ )
            {
              const double x = c[i];
              const double xm = cm[i];
              cm[i] = dt_inv * ( x - xm );
            }
            tmap["gram"].start();
            s_.wf.sd(ispin,ikp)->gram();
            tmap["gram"].stop();
          }
        }
      }
    }
    
    // compute ionic forces at last position for postprocessing of ionic
    // positions (Stoermer end step)
    double energy = 
      e.energy(false,dwf,compute_forces,fion,compute_stress,dcell); 
    ionic_stepper->postprocess(fion);
  }
  else
  {
    // delete wavefunction velocities
    if ( s_.wfv != 0 )
      delete s_.wfv;
    s_.wfv = 0;
  }

}
