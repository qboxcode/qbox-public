////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.5 2004-03-11 21:52:31 fgygi Exp $

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
#include "SDCellStepper.h"
#include "Preconditioner.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, EnergyFunctional& ef, int nite) : 
  SampleStepper(s), ef_(ef), 
  dwf(s.wf), wfv(s.wfv), nite_(nite) {}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(int niter)
{
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  const UnitCell& cell = wf.cell();
  const double omega = cell.volume();
  
  const double dt = s_.ctrl.dt;
  const double dt_inv = 1.0 / dt;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;
  
  const bool compute_hpsi = ( wf_dyn != "LOCKED" );
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  
  Timer tm_iter;
  
  const bool use_preconditioner = wf_dyn == "PSD" || wf_dyn == "PSDA";
  Preconditioner *preconditioner = 0;
  if ( use_preconditioner )
  {
    // create a preconditioner using the information about wf in s_.wf
    // and the information about the hessian in df
    preconditioner = new Preconditioner(s_,ef_);
  }
  
  WavefunctionStepper* wf_stepper = 0;
  if ( wf_dyn == "SD" )
    wf_stepper = new SDWavefunctionStepper(s_,tmap);
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(s_,*preconditioner,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(s_,*preconditioner,tmap);
  // wf_stepper == 0 indicates that wf_dyn == LOCKED
    
  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "MD" )
    ionic_stepper = new MDIonicStepper(s_);
    
  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);
  
  // Allocate wavefunction velocity if not available
  if ( atoms_dyn != "LOCKED" )
  {
    if ( s_.wfv == 0 )
    {
      s_.wfv = new Wavefunction(wf);
      s_.wfv->clear();
    }
  }
      
  if ( niter == 0 )
  {
    // evaluate and print energy
    double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
    if ( s_.ctxt_.onpe0() )
    {
      cout << ef_;
    }
  }
  
  for ( int iter = 0; iter < niter; iter++ )
  {
    // ionic iteration
 
    tm_iter.start();
 
    if ( s_.ctxt_.onpe0() )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
 
    double energy = 0.0;
    if ( compute_forces || compute_stress )
    {
      // compute energy and ionic forces using existing wavefunction
      energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

      if ( s_.ctxt_.onpe0() )
      {
        cout.setf(ios::fixed,ios::floatfield);
        cout.setf(ios::right,ios::adjustfield);
        cout << "  <ekin>   " << setprecision(8)
             << setw(15) << ef_.ekin() << " </ekin>\n";
        if ( compute_stress )
          cout << "  <econf>  " << setw(15) << ef_.econf() << " </econf>\n";
        cout << "  <eps>    " << setw(15) << ef_.eps() << " </eps>\n"
             << "  <enl>    " << setw(15) << ef_.enl() << " </enl>\n"
             << "  <ecoul>  " << setw(15) << ef_.ecoul() << " </ecoul>\n"
             << "  <exc>    " << setw(15) << ef_.exc() << " </exc>\n"
             << "  <esr>    " << setw(15) << ef_.esr() << " </esr>\n"
             << "  <eself>  " << setw(15) << ef_.eself() << " </eself>\n"
             << "  <etotal> " << setw(15) << ef_.etotal() << " </etotal>\n";
        if ( compute_stress )
        {
          const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
          const double enthalpy = ef_.etotal() + pext * cell.volume();
          cout << "  <enthalpy> " << setw(15) << enthalpy << " </enthalpy>\n"
             << flush;
        }
      }
    }
    
    // fion contains forces f0(r0)
 
    if ( compute_forces )
    {
      if ( iter > 0 ) 
      {
        ionic_stepper->compute_v0(fion);
        ionic_stepper->update_v();
      }
        
      ionic_stepper->compute_rp(fion);
      
      // Insert constraint enforcement here: using r0, rp, modify rp
      
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
                 << ionic_stepper->r0(is,i) << " "
                 << ionic_stepper->r0(is,i+1) << " " 
                 << ionic_stepper->r0(is,i+2) << " </position>\n"
                 << "    <velocity> " 
                 << ionic_stepper->v0(is,i) << " "
                 << ionic_stepper->v0(is,i+1) << " " 
                 << ionic_stepper->v0(is,i+2) << " </velocity>\n"
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
    }
    
    if ( compute_stress )
    {
      compute_sigma();
      print_stress();
      
      if ( cell_dyn != "LOCKED" )
      {
        cell_stepper->compute_new_cell(sigma);
 
        // Update cell
        cell_stepper->update_cell();
 
        ef_.cell_moved();
        ef_.atoms_moved(); // modifications of the cell also move ions
        
        if ( use_preconditioner )
          preconditioner->update();
      }
    }
    
    if ( compute_forces )
    {
      ionic_stepper->update_r();
      // velocities in ionic_stepper are v(t-dt)
      ef_.atoms_moved();
    }
    
    // Recalculate ground state wavefunctions
    
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
        double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
 
        wf_stepper->update(dwf);
 
        if ( s_.ctxt_.onpe0() )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <ekin_int>   " << setprecision(8)
               << setw(15) << ef_.ekin() << " </ekin_int>\n";
          if ( compute_stress )
          {
            cout << "  <econf_int>  " << setw(15) << ef_.econf() 
                 << " </econf_int>\n";
          }
          cout << "  <eps_int>    " << setw(15) << ef_.eps() << " </eps_int>\n"
               << "  <enl_int>    " << setw(15) << ef_.enl() << " </enl_int>\n"
               << "  <ecoul_int>  " << setw(15) << ef_.ecoul() << " </ecoul_int>\n"
               << "  <exc_int>    " << setw(15) << ef_.exc() << " </exc_int>\n"
               << "  <esr_int>    " << setw(15) << ef_.esr() << " </esr_int>\n"
               << "  <eself_int>  " << setw(15) << ef_.eself() << " </eself_int>\n"
               << "  <etotal_int> " << setw(15) << ef_.etotal() << " </etotal_int>\n";
          if ( compute_stress )
          {
            const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
            const double enthalpy = ef_.etotal() + pext * cell.volume();
            cout << "  <pv> " << setw(15) << pext * cell.volume()
                 << " </pv>" << endl;
            cout << "  <enthalpy_int> " << setw(15) << enthalpy << " </enthalpy_int>\n"
                 << flush;
          }
        }
      }
      wf_stepper->postprocess();
    }
    else
    {
      // wf_stepper == 0, wf_dyn == LOCKED
      // evaluate and print energy
      double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
      if ( s_.ctxt_.onpe0() )
      {
        cout << ef_;
      }
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
    
    // compute ionic forces at last position to update velocities
    // consistently with last position
    double energy = 
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
      
    // Note: next line function call updates velocities in the AtomSet
    ionic_stepper->compute_v0(fion);
    ionic_stepper->update_v();
  }
  else
  {
    // delete wavefunction velocities
    if ( s_.wfv != 0 )
      delete s_.wfv;
    s_.wfv = 0;
  }
  
  // delete steppers
  if ( wf_stepper != 0 ) delete wf_stepper;
  if ( ionic_stepper != 0 ) delete ionic_stepper;
  if ( cell_stepper != 0 ) delete cell_stepper;
  
  // delete preconditioner
  if ( use_preconditioner ) delete preconditioner;
}
