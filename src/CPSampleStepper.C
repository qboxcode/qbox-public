////////////////////////////////////////////////////////////////////////////////
//
// CPSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.C,v 1.8 2004-10-28 16:54:28 fgygi Exp $

#include "CPSampleStepper.h"
#include "SlaterDet.h"
#include "MDWavefunctionStepper.h"
#include "MDIonicStepper.h"
#include "SDCellStepper.h"
#include "Basis.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CPSampleStepper::CPSampleStepper(Sample& s) : 
  SampleStepper(s), cd_(s.wf), ef_(s,cd_), dwf(s.wf), wfv(s.wfv)
{
  mdwf_stepper = new MDWavefunctionStepper(s,tmap);
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
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void CPSampleStepper::step(int niter)
{
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  
  const double dt = s_.ctrl.dt;
  double ekin_ion=0.0,ekin_e, temp_ion, eta;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  assert(wf_dyn=="MD");
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;
  
  const bool compute_hpsi = true;
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);
    
  if ( s_.wfv == 0 )
  {
    s_.wfv = new Wavefunction(wf);
    s_.wfv->clear();
  }
  
  Timer tm_iter;
  
  tmap["charge"].start();
  cd_.update_density();
  tmap["charge"].stop();
  
  ef_.update_vhxc();
  double energy =
    ef_.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma_eks);
 
  mdwf_stepper->compute_wfm(dwf);

  for ( int iter = 0; iter < niter; iter++ )
  {
    tm_iter.reset();
    tm_iter.start();
    if ( s_.ctxt_.mype() == 0 )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
 
    mdwf_stepper->update(dwf);
      
    ekin_e = mdwf_stepper->ekin();
 
    if ( s_.ctxt_.onpe0() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout.setf(ios::right,ios::adjustfield);
      cout << "  <ekin>     " << setprecision(8)
           << setw(15) << ef_.ekin() << " </ekin>\n";
      if ( use_confinement )
      {
        cout << "  <econf>    " << setw(15) << ef_.econf()
             << " </econf>\n";
      }
      cout << "  <eps>      " << setw(15) << ef_.eps() << " </eps>\n"
           << "  <enl>      " << setw(15) << ef_.enl() << " </enl>\n"
           << "  <ecoul>    " << setw(15) << ef_.ecoul() << " </ecoul>\n"
           << "  <exc>      " << setw(15) << ef_.exc() << " </exc>\n"
           << "  <esr>      " << setw(15) << ef_.esr() << " </esr>\n"
           << "  <eself>    " << setw(15) << ef_.eself() << " </eself>\n"
           << "  <etotal>   " << setw(15) << ef_.etotal() << " </etotal>\n"
           << flush;
      if ( compute_stress )
      {
        const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
        const double enthalpy = ef_.etotal() + pext * s_.wf.cell().volume();
        cout << "  <pv>     " << setw(15) << pext * s_.wf.cell().volume()
             << " </pv>" << endl;
        cout << "  <enthalpy> " << setw(15) << enthalpy << " </enthalpy>\n"
           << flush;
      }
    }
 
    if ( compute_forces )
    {
      if ( iter > 0 ) 
      {
        mdionic_stepper->compute_v0(fion);
        mdionic_stepper->update_v();
      }
      mdionic_stepper->compute_rp(fion);
      ekin_ion = mdionic_stepper->ekin();
      
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
      }
    }
 
    if ( s_.ctxt_.onpe0() )
    {
      cout << "  <ekin_e> " << ekin_e << " </ekin_e>\n";
 
      if ( compute_forces )
      {
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
        cout << "  <temp_ion> " << mdionic_stepper->temp() << " </temp_ion>\n";
        cout << "  <eta_ion> " << mdionic_stepper->eta() << " </eta_ion>\n";
      }
      cout << "  <econst> " << energy+ekin_ion+ekin_e << " </econst>\n";
      cout << "  <ekin_ec> " << energy+ekin_ion+2*ekin_e << " </ekin_ec>\n";
    }
    
    if ( compute_stress )
    {
      if ( s_.ctxt_.onpe0() )            
      {                                  
        cout << "<unit_cell>" << endl;   
        cout << s_.wf.cell();            
        cout << "</unit_cell>" << endl;  
      }                                  
      compute_sigma();
      print_stress();
      
      if ( cell_dyn != "LOCKED" )
      {
        cell_stepper->compute_new_cell(sigma);
      
        // Update cell
        cell_stepper->update_cell();
      
        ef_.cell_moved();
        ef_.atoms_moved();
      }
    }
    
    if ( compute_forces )
    {
      mdionic_stepper->update_r();
      ef_.atoms_moved();
    }
    
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    ef_.update_vhxc();
    energy =
      ef_.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma_eks);

    if ( s_.ctxt_.mype() == 0 )
      cout << "  </iteration>" << endl;

    // print iteration time
    tm_iter.stop();
    double time = tm_iter.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "  <!-- timing "
           << setw(15) << "iteration"
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  } // iter
 
  // dwf and fion now contain the forces on wavefunctions and ions at the
  // endpoint
 
  mdwf_stepper->compute_wfv(dwf); // replace wfm by wfv
  
  if ( compute_forces )
  {
    // Note: next line function call updates velocities in the AtomSet
    mdionic_stepper->compute_v0(fion);
    mdionic_stepper->update_v();
  }
  
  if ( cell_stepper != 0 ) delete cell_stepper;
}
