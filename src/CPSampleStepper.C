////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: CPSampleStepper.C,v 1.2 2004-02-04 19:55:16 fgygi Exp $

#include "CPSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "MDWavefunctionStepper.h"
#include "MDIonicStepper.h"
#include "Basis.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CPSampleStepper::CPSampleStepper(Sample& s) : SampleStepper(s), dwf(s.wf), 
  wfv(s.wfv)
{
  mdwf_stepper = new MDWavefunctionStepper(s,tmap);
  assert(mdwf_stepper!=0);
  mdionic_stepper = 0;
  if ( s.ctrl.atoms_dyn != "LOCKED" )
    mdionic_stepper = new MDIonicStepper(s);
  fion.resize(s_.atoms.nsp());
  for ( int is = 0; is < fion.size(); is++ )
    fion[is].resize(3*s_.atoms.na(is));
}

////////////////////////////////////////////////////////////////////////////////
void CPSampleStepper::step(EnergyFunctional& e, int niter)
{
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  valarray<double> sigma(6);
  
  const double dt = s_.ctrl.dt;
  double ekin_ion=0.0,ekin_e, temp_ion, eta;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  assert(wf_dyn=="MD");
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  
  const bool compute_hpsi = true;
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = false;

  Timer tm_iter;
  
  double energy =
    e.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma);
 
  for ( int iter = 0; iter < niter; iter++ )
  {
    tm_iter.reset();
    tm_iter.start();
    if ( s_.ctxt_.mype() == 0 )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
 
    if ( iter == 0 )
      mdwf_stepper->stoermer_start(dwf);
    else
      mdwf_stepper->update(dwf);
      
    ekin_e = mdwf_stepper->ekin();
 
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
                 << mdionic_stepper->tau0(is,i) << " "
                 << mdionic_stepper->tau0(is,i+1) << " " 
                 << mdionic_stepper->tau0(is,i+2) << " </position>\n"
                 << "    <force> " << fion[is][i] << " "
                 << fion[is][i+1] << " " << fion[is][i+2]
                 << " </force>\n  </atom>" << endl;
 
            i += 3;
          }
        }
      }
      
      if ( iter == 0 )
        mdionic_stepper->preprocess(fion);
        
      mdionic_stepper->update(fion);
        
      ekin_ion = mdionic_stepper->ekin();
      e.atoms_moved();
    }
 
    if ( s_.wf.context().onpe0() )
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
    
    energy =
      e.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma);

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
 
  mdwf_stepper->stoermer_end(dwf); // replace wfm by wfv
  
  if ( compute_forces )
  {
    mdionic_stepper->postprocess(fion);
  }
}
