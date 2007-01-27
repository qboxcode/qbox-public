////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.30 2007-01-27 23:46:30 fgygi Exp $

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
#include "SDIonicStepper.h"
#include "SDAIonicStepper.h"
#include "MDIonicStepper.h"
#include "SDCellStepper.h"
#include "Preconditioner.h"
#include "AndersonMixer.h"

#ifdef USE_APC
#include "apc.h"
#endif

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nitscf, int nite) : 
  SampleStepper(s), cd_(s.wf), ef_(s,cd_), 
  dwf(s.wf), wfv(s.wfv), nitscf_(nitscf), nite_(nite) {}

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::~BOSampleStepper()
{
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
void BOSampleStepper::step(int niter)
{
  const bool onpe0 = s_.ctxt_.onpe0();
  // determine whether eigenvectors must be computed
  // eigenvectors are computed if explicitly requested with wf_diag==T
  // or if the SlaterDet has fractionally occupied states
  const bool fractional_occ = (s_.wf.nel() != 2 * s_.wf.nst());
  const bool compute_eigvec = fractional_occ || s_.ctrl.wf_diag == "T";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };
  
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  const int nspin = wf.nspin();
  
  const UnitCell& cell = wf.cell();
  const double omega = cell.volume();
  
  const double dt = s_.ctrl.dt;
  const double dt_inv = 1.0 / dt;
  
  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;
  
  const bool extrapolate_wf = true;
  const bool ntc_extrapolation = 
    s_.ctrl.debug.find("NTC_EXTRAPOLATION") != string::npos;
  const bool asp_extrapolation = 
    s_.ctrl.debug.find("ASP_EXTRAPOLATION") != string::npos;

  Wavefunction* wfmm;
  if ( ntc_extrapolation || asp_extrapolation ) wfmm = new Wavefunction(wf);
  
  const bool compute_hpsi = ( wf_dyn != "LOCKED" );
  const bool compute_forces = ( atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool use_confinement = ( s_.ctrl.ecuts > 0.0 );
  
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
  {
    const double emass = s_.ctrl.emass;
    double dt2bye = (emass == 0.0) ? 0.5 / wf.ecut() : dt*dt/emass;
  
    // divide dt2bye by facs coefficient if stress == ON
    const double facs = 2.0;
    if ( s_.ctrl.stress == "ON" )
    {
      dt2bye /= facs;
    }
    wf_stepper = new SDWavefunctionStepper(wf,dt2bye,tmap);
  }
  else if ( wf_dyn == "PSD" )
    wf_stepper = new PSDWavefunctionStepper(wf,*preconditioner,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(wf,*preconditioner,tmap);
    
  // wf_stepper == 0 indicates that wf_dyn == LOCKED
    
  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "SDA" )
    ionic_stepper = new SDAIonicStepper(s_);
  else if ( atoms_dyn == "MD" )
    ionic_stepper = new MDIonicStepper(s_);
    
  if ( ionic_stepper )
    ionic_stepper->setup_constraints();
    
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
    const bool compute_hpsi = false;
    const bool compute_forces = false;
    const bool compute_stress = false;
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    ef_.update_vhxc();
    double energy = ef_.energy(compute_hpsi,dwf,compute_forces,fion,
                               compute_stress,sigma_eks);
    if ( onpe0 )
    {
      cout << ef_;
    }
  }
  
  for ( int iter = 0; iter < niter; iter++ )
  {
    // ionic iteration
 
    tm_iter.start();
#ifdef USE_APC
    ApcStart(1);
#endif
 
    if ( onpe0 )
      cout << "  <iteration count=\"" << iter+1 << "\">\n";
      
    if ( ionic_stepper )
      atoms.sync();
 
    double energy = 0.0;
    if ( compute_forces || compute_stress )
    {
      // compute energy and ionic forces using existing wavefunction

      tmap["charge"].start();
      cd_.update_density();
      tmap["charge"].stop();
      
      ef_.update_vhxc();
      energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

      if ( onpe0 )
      {
        cout.setf(ios::fixed,ios::floatfield);
        cout.setf(ios::right,ios::adjustfield);
        cout << "  <ekin>   " << setprecision(8)
             << setw(15) << ef_.ekin() << " </ekin>\n";
        if ( use_confinement )
          cout << "  <econf>  " << setw(15) << ef_.econf() << " </econf>\n";
        cout << "  <eps>    " << setw(15) << ef_.eps() << " </eps>\n"
             << "  <enl>    " << setw(15) << ef_.enl() << " </enl>\n"
             << "  <ecoul>  " << setw(15) << ef_.ecoul() << " </ecoul>\n"
             << "  <exc>    " << setw(15) << ef_.exc() << " </exc>\n"
             << "  <esr>    " << setw(15) << ef_.esr() << " </esr>\n"
             << "  <eself>  " << setw(15) << ef_.eself() << " </eself>\n"
             << "  <ets>    " << setw(15) << ef_.ets() << " </ets>\n"
             << "  <etotal> " << setw(15) << ef_.etotal() << " </etotal>\n";
        if ( compute_stress )
        {
          const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
          const double enthalpy = ef_.etotal() + pext * cell.volume();
          cout << "  <pv>     " << setw(15) << pext * cell.volume()
               << " </pv>" << endl;
          cout << "  <enthalpy> " << setw(15) << enthalpy << " </enthalpy>\n"
             << flush;
        }
      }
    }
    
    if ( compute_forces )
    {
      if ( iter > 0 ) 
      {
        ionic_stepper->compute_v(energy,fion);
      }
      // at this point, positions r0, velocities v0 and forces fion are
      // consistent
      const double ekin_ion = ionic_stepper->ekin();                 
      const double temp_ion = ionic_stepper->temp();                 
 
      // print positions, velocities and forces at time t0           
      if ( onpe0 )
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
                 << " </force>\n";                                 
            cout << "  </atom>" << endl;                             
            i += 3;                                                  
          }                                                          
        }
        cout << "  <econst> " << energy+ekin_ion << " </econst>\n";  
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";     
        cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";     
      }
    }
    
    if ( compute_forces )
    {
      if ( s_.constraints.size() > 0 )
      {
        s_.constraints.compute_forces(ionic_stepper->r0(), fion);
        if ( onpe0 )
        {
          s_.constraints.list_constraints(cout);
        }
      }
      // move atoms to new position: r0 <- r0 + v0*dt + dt2/m * fion
      ionic_stepper->compute_r(energy,fion);
      ef_.atoms_moved();
    }
    
    if ( compute_stress )
    {
      if ( onpe0 )            
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
        ef_.atoms_moved(); // modifications of the cell also move ions
        
        if ( use_preconditioner )
          preconditioner->update();
      }
    }
    
    // Recalculate ground state wavefunctions
    
    // wavefunction extrapolation
    if ( compute_forces && extrapolate_wf )
    {
      for ( int ispin = 0; ispin < nspin; ispin++ )
      {
        for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
        {
          if ( s_.wf.sd(ispin,ikp) != 0 )
          {
            if ( s_.wf.sdcontext(ispin,ikp)->active() )
            {
              if ( ntc_extrapolation )
              {
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();      
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();   
                double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();     
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();              
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                const int len = 2*mloc*nloc;              
                if ( iter == 0 )
                {
                  for ( int i = 0; i < len; i++ )
                  {
                    const double x = c[i];
                    const double v = cv[i];
                    // extrapolation using velocity in cv
                    c[i] = x + dt * v;
                    cv[i] = x;                                        
                  }
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                }
                else if ( iter == 1 )
                {
                  s_.wfv->align(s_.wf);
                  for ( int i = 0; i < len; i++ )
                  {
                    const double x = c[i];
                    const double xm = cv[i];
                    c[i] = 2.0 * x - xm;
                    cv[i] = x;                             
                    cmm[i] = xm;
                  }
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                }
                else
                {
                  // align wf with wfmm before extrapolation
                  // s_.wf.align(*wfmm);
                  wfmm->align(s_.wf);
                
                  // extrapolate
                  for ( int i = 0; i < len; i++ )                        
                  {                                                            
                    const double x = c[i];   // current wf (scf converged) at t
                    const double xm = cv[i]; // extrapolated wf at t           
                    const double xmm = cmm[i]; // extrapolated wf at t-dt 
                    c[i] = 2.0 * x - xmm;                                      
                    // save extrapolated value at t in cmm                             
                    cmm[i] = xm;           
                  }
                  
                  // orthogonalize the extrapolated value        
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                  //tmap["lowdin"].start();                                    
                  //s_.wf.sd(ispin,ikp)->lowdin();                             
                  //tmap["lowdin"].stop();     
                                                  
                  // c[i] now contains the extrapolated value      
                  // save a copy in cv[i]                          
                  for ( int i = 0; i < len; i++ )                        
                  {                                                            
                    cv[i] = c[i];                                        
                  }                                                
                }
                // c[i] is now ready for electronic iterations
              }
              else if ( asp_extrapolation )
              {
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();      
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();   
                double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();     
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();              
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                const int len = 2*mloc*nloc;              
                if ( iter == 0 )
                {
                  for ( int i = 0; i < len; i++ )
                  {
                    const double x = c[i];
                    const double v = cv[i];
                    // extrapolation using velocity in cv
                    c[i] = x + dt * v;
                    cv[i] = x;                                        
                  }
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                }
                else if ( iter == 1 )
                {
                  //s_.wfv->align(s_.wf);
                  for ( int i = 0; i < len; i++ )
                  {
                    const double x = c[i];
                    const double xm = cv[i];
                    c[i] = 2.0 * x - xm;
                    cv[i] = x;                             
                    cmm[i] = xm;
                  }
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                }
                else
                {
                  // align wf with wfmm before extrapolation
                  // s_.wf.align(*wfmm);
                  // wfmm->align(s_.wf);
                
                  // extrapolate
                  for ( int i = 0; i < len; i++ )                        
                  {                                                            
                    const double x = c[i];   // current wf (scf converged) at t
                    const double xm = cv[i]; // extrapolated wf at t           
                    const double xmm = cmm[i]; // extrapolated wf at t-dt 
                    const double asp_a1 = 0.5;
                    c[i] = 2.0 * x - xm + 
                           asp_a1 * ( x - 2.0 * xm + xmm );                                      
                    //c[i] = 2.5 * x - 2.0 * xm + 0.5 * xmm;                                      
                    cmm[i] = xm;
                    cv[i] = x;           
                  }
                  
                  // orthogonalize the extrapolated value        
                  tmap["gram"].start();                          
                  s_.wf.sd(ispin,ikp)->gram();                   
                  tmap["gram"].stop();                           
                  //tmap["lowdin"].start();                                    
                  //s_.wf.sd(ispin,ikp)->lowdin();                             
                  //tmap["lowdin"].stop();     
                                                  
                  // c[i] now contains the extrapolated value      
                }
                // c[i] is now ready for electronic iterations
              }
              else // normal extrapolation
              {
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();      
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();   
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();              
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                const int len = 2*mloc*nloc;              
                if ( iter == 0 )
                {
                  // replace velocity in cv by cm = c - dt * cv
                  for ( int i = 0; i < len; i++ )
                  {
                    const double x = c[i];
                    const double v = cv[i];
                    cv[i] = x - dt * v;
                  }
                }
                // linear extrapolation                                       
                for ( int i = 0; i < len; i++ )                        
                {                                                              
                  const double x = c[i];
                  const double xm = cv[i];
                  c[i] = 2.0 * x - xm;                                       
                  cv[i] = x;                                                 
                }                                                              
                tmap["gram"].start();                                        
                s_.wf.sd(ispin,ikp)->gram();                                 
                tmap["gram"].stop();                                         
 
                //tmap["lowdin"].start();                                      
                //s_.wf.sd(ispin,ikp)->lowdin();                               
                //tmap["lowdin"].stop();                                       
 
                //tmap["ortho_align"].start();                                 
                //s_.wf.sd(ispin,ikp)->ortho_align(*s_.wfv->sd(ispin,ikp));    
                //tmap["ortho_align"].stop();                                  
 
                //tmap["riccati"].start();                                     
                //s_.wf.sd(ispin,ikp)->riccati(*s_.wfv->sd(ispin,ikp));        
                //tmap["riccati"].stop();                                      
              }
            } // if active
          }
        }
      }
    } // compute_forces && extrapolate_wf

    // do nitscf self-consistent iterations, each with nite electronic steps
    if ( wf_stepper != 0 )
    {
      assert(cd_.rhog.size()==1); // works for nspin=1 only
      vector<complex<double> > rhog_current(cd_.rhog[0]);
      vector<complex<double> > rhog_last(rhog_current);
      vector<complex<double> > drhog(rhog_current.size());
      vector<complex<double> > drhog_bar(rhog_current.size());
      AndersonMixer mixer(2*rhog_current.size(),&cd_.vcontext());
      mixer.set_theta_max(2.0);
      mixer.set_theta_nc(1.0);
      
      wf_stepper->preprocess();
      for ( int itscf = 0; itscf < nitscf_; itscf++ )
      {
        if ( nite_ > 1 && onpe0 )
          cout << "  <!-- BOSampleStepper: start scf iteration -->" << endl;
          
        // compute new density in cd_.rhog
        tmap["charge"].start();
        cd_.update_density();
        tmap["charge"].stop();
        
        // charge mixing
        if ( nite_ > 1 )
        {
          if ( itscf == 0 )
          {
            rhog_current = cd_.rhog[0];
          }

          // compute correction drhog
          for ( int i=0; i < rhog_current.size(); i++ )
          {
            drhog[i] = (cd_.rhog[0][i] - rhog_current[i]);
          }

          // Apply Kerker preconditioner to drhog                        
          // Use Kerker preconditioning if rc_Kerker > 0.0,                  
          // no preconditioning otherwise                                    
          const double alpha = s_.ctrl.charge_mix_coeff;                     
          const double *const g2 = cd_.vbasis()->g2_ptr();                   
          // real space Kerker cutoff in a.u.                                
          const double rc_Kerker = s_.ctrl.charge_mix_rcut;                  
          if ( rc_Kerker > 0.0 )                                             
          {                                                                  
            const double q0_kerker = 2 * M_PI / rc_Kerker;                   
            const double q0_kerker2 = q0_kerker * q0_kerker;                 
            for ( int i=0; i < rhog_current.size(); i++ )                    
            {                                                                
              drhog[i] *= alpha * g2[i] / ( g2[i] + q0_kerker2 );        
            }                                                                
          }                                                                  
          else                                                               
          {                                                                  
            for ( int i=0; i < rhog_current.size(); i++ )                    
            {                                                                
              drhog[i] *= alpha;                                         
            }                                                                
          }

          // Anderson acceleration
          double theta = 0.0;
          for ( int i=0; i < rhog_current.size(); i++ )  
          {                                              
            drhog_bar[i] = drhog[i];                     
          }                                              

          const bool anderson_charge_mixing = true;
          if ( anderson_charge_mixing )
          {
            mixer.update((double*)&drhog[0],&theta,(double*)&drhog_bar[0]);
            if ( onpe0 )
            {
              cout << "  <!-- Charge mixing: Anderson theta="
                   << theta << " -->" << endl;
            }
          }
          
          // update rhog_current 
          // rhog_current = rhog_current + theta*(rhog_current-rhog_last)
 
          for ( int i=0; i < rhog_current.size(); i++ )
          {
            complex<double> rhotmp = rhog_current[i];
            rhog_current[i] += theta * (rhog_current[i] - rhog_last[i]);
            rhog_last[i] = rhotmp;
          }

          // Apply correction
          for ( int i=0; i < rhog_current.size(); i++ )
          {
            cd_.rhog[0][i] = rhog_current[i] + drhog_bar[i]; 
          }
          rhog_current = cd_.rhog[0];
          cd_.update_rhor();
        } // if nite > 1

        ef_.update_vhxc();
        
        // reset stepper only if multiple non-selfconsistent steps
        if ( nite_ > 1 ) wf_stepper->preprocess();
      
        for ( int ite = 0; ite < nite_; ite++ )
        {
          double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
          
          // compute the sum of eigenvalues (with fixed weight)
          // to measure convergence of the subspace update
          // compute trace of the Hamiltonian matrix Y^T H Y
          // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
          const double eigenvalue_sum = s_.wf.dot(dwf);
          if ( onpe0 )
            cout << "  <eigenvalue_sum> "
                 << eigenvalue_sum << " </eigenvalue_sum>" << endl;
 
          wf_stepper->update(dwf);
          
          if ( onpe0 )
          {
            cout.setf(ios::fixed,ios::floatfield);
            cout.setf(ios::right,ios::adjustfield);
            cout << "  <etotal_int> " << setw(15)
                 << energy << " </etotal_int>\n";
            if ( compute_stress )
            {
              const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
              const double enthalpy = energy + pext * cell.volume();
              cout << "  <enthalpy_int> " << setw(15) 
                   << enthalpy << " </enthalpy_int>\n"
                   << flush;
            }
          }
        } // for ite
        
        // subspace diagonalization
        if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" )
        {
          energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
          s_.wf.diag(dwf,compute_eigvec);
          if ( onpe0 )                                    
          {
            // print eigenvalues
            for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
            {
              for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
              {
                if ( wf.sd(ispin,ikp) != 0 )
                {
                  if ( wf.sdcontext(ispin,ikp)->active() )
                  {
                    const int nst = wf.sd(ispin,ikp)->nst();
                    const double eVolt = 2.0 * 13.6058;                 
                    cout <<    "  <eigenvalues spin=\"" << ispin        
                         << "\" kpoint=\"" << wf.sd(ispin,ikp)->kpoint()   
                         << "\" n=\"" << nst << "\">" << endl;          
                    for ( int i = 0; i < nst; i++ )                   
                    {                                                   
                      cout << setw(12) << setprecision(5) 
                           << wf.sd(ispin,ikp)->eig(i)*eVolt;
                      if ( i%5 == 4 ) cout << endl;                     
                    }                                                   
                    if ( nst%5 != 0 ) cout << endl;                   
                    cout << "  </eigenvalues>" << endl;
                  }
                }
              }
            }                 
          }                                                       
        }
        
        // update occupation numbers if fractionally occupied states
        if ( fractional_occ )
        {
          wf.update_occ(s_.ctrl.fermi_temp);
          const double wf_entropy = wf.entropy();
          if ( onpe0 )
          {
            cout << "  <!-- Wavefunction entropy: " << wf_entropy
                 << " -->" << endl;
            const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
            cout << "  <!-- Entropy contribution to free energy: "
                 << - wf_entropy * s_.ctrl.fermi_temp * boltz
                 << " -->" << endl;
          }
        }
        
        if ( nite_ > 1 && onpe0 )
          cout << "  <!-- BOSampleStepper: end scf iteration -->" << endl;
      } // for itscf
      
      
      if ( !compute_forces && !compute_stress )
        if ( onpe0 )
          cout << ef_;
          
      wf_stepper->postprocess();
    }
    else
    {
      // wf_stepper == 0, wf_dyn == LOCKED
      // evaluate and print energy
      tmap["charge"].start();
      cd_.update_density();
      tmap["charge"].stop();
      ef_.update_vhxc();
      double energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
      if ( onpe0 )
      {
        cout << ef_;
      }
    }
 
#ifdef USE_APC
    ApcStop(1);
#endif
    // print iteration time
    double time = tm_iter.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( onpe0 )
    {
      cout << "  <!-- timing "
           << setw(15) << "iteration"
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
      cout << "  </iteration>" << endl;
    }
    if ( compute_forces )
      s_.constraints.update_constraints(dt);
  } // for iter
  
  if ( compute_forces )
  {
    // compute wavefunction velocity after last iteration
    // s_.wfv contains the previous wavefunction
    
    // if eigenvectors were computed, use alignment before computing velocity
    if ( compute_eigvec )
    {
      s_.wfv->align(s_.wf);
    }
      
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
            const int len = 2*mloc*nloc;              
            if ( ntc_extrapolation )
            {
              double* cmm = (double*) wfmm->sd(ispin,ikp)->c().cvalptr();     
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xmm = cmm[i];
                cm[i] = dt_inv * ( x - xmm );
              }
              tmap["gram"].start();
              s_.wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
            else // normal extrapolation or asp_extrapolation
            {
              for ( int i = 0; i < len; i++ )
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
    }
    
    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();
    
    ef_.update_vhxc();
    double energy = 
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
      
    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
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
  if ( ntc_extrapolation || asp_extrapolation ) delete wfmm;
}
