////////////////////////////////////////////////////////////////////////////////
//
// BOSampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.C,v 1.12 2004-09-15 01:01:28 fgygi Exp $

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

#define POTENTIAL_MIXING 1
#define CHARGE_MIXING 0

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nitscf, int nite) : 
  SampleStepper(s), cd_(s.wf), ef_(s,cd_), 
  dwf(s.wf), wfv(s.wfv), nitscf_(nitscf), nite_(nite) {}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(int niter)
{
  const bool extrapolate_wf = true;
  const bool quad_extrapolation = false;
  const bool compute_fv = false;
  const int nempty = s_.wf.nempty();
  const bool compute_eigvec = nempty > 0 || s_.ctrl.wf_diag == "T";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };
  
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  Wavefunction wfmm(wf);
  const int nspin = wf.nspin();
  
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
    const bool compute_hpsi = false;
    const bool compute_forces = false;
    const bool compute_stress = false;
    cd_.update_density();
    ef_.update_vhxc();
    double energy = ef_.energy(compute_hpsi,dwf,compute_forces,fion,
                               compute_stress,sigma_eks);
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
    
      cd_.update_density();
      ef_.update_vhxc();
      energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);

      if ( s_.ctxt_.onpe0() )
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
             << "  <etotal> " << setw(15) << ef_.etotal() << " </etotal>\n";
        if ( compute_stress )
        {
          const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
          const double enthalpy = ef_.etotal() + pext * cell.volume();
          cout << "  <pv> " << setw(15) << pext * cell.volume()
               << " </pv>" << endl;
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
    if ( compute_forces && extrapolate_wf )
    {
      // extrapolate wavefunctions
      // s_.wfv contains the wavefunction velocity
      
      if ( iter > 0 && compute_eigvec )
      {
        // eigenvectors were computed, need alignment
        // wfv contains wfm since iter > 0
        s_.wfv->align(s_.wf);
      }
      
      for ( int ispin = 0; ispin < nspin; ispin++ )
      {
        for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
        {
          if ( s_.wf.sd(ispin,ikp) != 0 )
          {
            if ( s_.wf.sdcontext(ispin,ikp)->active() )
            {
              ortho_type ortho;
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
                //ortho = GRAM; 
                //ortho = ORTHO_ALIGN;
                ortho = LOWDIN;
              }
              else if ( iter == 1  || !quad_extrapolation )
              {
                // linear extrapolation using the previous wavefunction
                // cm is stored in wfv
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                double* cmm = (double*) wfmm.sd(ispin,ikp)->c().cvalptr();
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                for ( int i = 0; i < 2*mloc*nloc; i++ )
                {
                  const double x = c[i];
                  const double xm = cv[i];
                  c[i] = 2.0 * x - xm;
                  cv[i] = x;
                  cmm[i] = xm;
                }
                ortho = ORTHO_ALIGN;
                //ortho = RICCATI;
                //ortho = LOWDIN;
                //ortho = GRAM;
              }
              else
              {
                // quad_extrapolation == true && iter >= 1
                // quadratic extrapolation using the two previous wavefunctions
                // cm is stored in wfv, cmm is stored in wfmm
                double* c = (double*) s_.wf.sd(ispin,ikp)->c().cvalptr();
                double* cv = (double*) s_.wfv->sd(ispin,ikp)->c().cvalptr();
                double* cmm = (double*) wfmm.sd(ispin,ikp)->c().cvalptr();
                const int mloc = s_.wf.sd(ispin,ikp)->c().mloc();
                const int nloc = s_.wf.sd(ispin,ikp)->c().nloc();
                for ( int i = 0; i < 2*mloc*nloc; i++ )
                {
                  const double x = c[i];
                  const double xm = cv[i];
                  const double xmm = cmm[i];
                  
                  // linear extrapolation
                  //c[i] = 2.0 * x - xm;
                  
                  // Lagrange quadratic extrapolation
                  //c[i] = 3.0 * ( x - xm ) + xmm;
                  
                  // k=1 extrapolation 
                  c[i] = 2.0 * x - xm + 0.5 * ( x - 2.0 * xm + xmm );
                  
                  // damped/enhanced linear extrapolation
                  //c[i] = x + 0.9 * ( x - xm );
                  
                  cv[i] = x;
                  cmm[i] = xm;
                }
                ortho = ORTHO_ALIGN;
                //ortho = RICCATI;
                //ortho = LOWDIN;
                //ortho = GRAM;
              }
              
              switch ( ortho )
              {
                case GRAM:
                  tmap["gram"].stop();
                  s_.wf.sd(ispin,ikp)->gram();
                  tmap["gram"].stop();
                  break;
 
                case LOWDIN:
                  tmap["lowdin"].stop();
                  s_.wf.sd(ispin,ikp)->lowdin();
                  tmap["lowdin"].stop();
                  break;
 
                case ORTHO_ALIGN:
                  tmap["ortho_align"].stop();
                  s_.wf.sd(ispin,ikp)->ortho_align(*s_.wfv->sd(ispin,ikp));
                  tmap["ortho_align"].stop();
                  break;
 
                case RICCATI:
                  tmap["riccati"].stop();
                  s_.wf.sd(ispin,ikp)->riccati(*s_.wfv->sd(ispin,ikp));
                  tmap["riccati"].stop();
                  break;
              }
            } // if active
          }
        }
      }
    } // compute_forces && extrapolate_wf

    // do nitscf self-consistent iterations, each with nite electronic steps
    if ( wf_stepper != 0 )
    {
#if CHARGE_MIXING
      vector<vector<complex<double> > > rhog_old(cd_.rhog);
#endif
#if POTENTIAL_MIXING
      vector<vector<double> > vi(ef_.v_r);
      vector<vector<double> > vo1(vi), vi1(vi);
#endif
      
      for ( int itscf = 0; itscf < nitscf_; itscf++ )
      {
        if ( s_.ctxt_.onpe0() )
          cout << "  <!-- BOSampleStepper: start scf iteration -->" << endl;
        cd_.update_density();
        
        // charge mixing
        
        // compute correction delta_rhog = rhog_new - rhog_old
        // precondition with Kerker preconditioner
        
#if CHARGE_MIXING
        if ( itscf > 0 )
        {
          // mix density with previous density rhog_old
          const double alpha = 0.5;
          // Kerker cutoff: 15 a.u. in real space
          const double q0_kerker = 2.0 * M_PI / 15.0;
          const double q0_kerker2 = q0_kerker * q0_kerker;
          assert(rhog_old.size()==1); // assume nspin==1 for now
          const double *const g2 = cd_.vbasis()->g2_ptr();
          for ( int i=0; i < rhog_old[0].size(); i++ )
          {
            const complex<double> drhog = cd_.rhog[0][i] - rhog_old[0][i];
            const double fac = g2[i] / ( g2[i] + q0_kerker2 );
            cd_.rhog[0][i] = rhog_old[0][i] + alpha * fac * drhog;
          }
          cd_.update_rhor();
        }
        rhog_old = cd_.rhog;
#endif

        ef_.update_vhxc();
        
        // potential mixing
        // vlocal_old, vlocal_new in ef_
        
#if POTENTIAL_MIXING
        if ( nite_ > 1 )
        {
        
        // Potential mixing using Hamann's implementation of Anderson's method
        // generate next iteration using d. g. anderson's method
        assert(nspin==1);
        vector<vector<double> >& vo = ef_.v_r;
        const int size = vo[0].size();
        double thl = 0.0;
        if ( itscf == 0 )
        {
          // first iteration: only vo is defined. Use vo.
          for ( int i = 0; i < size; i++ )
          {
            const double bl = 0.3;
            double vn = vo[0][i];
            vi1[0][i] = vi[0][i];
            vo1[0][i] = vo[0][i];
            vi[0][i]  = vn;
            
          }
        }
        else if ( itscf == 1 )
        {
          // second iteration: only vo, vi are defined. Use thl = 0.0
          for ( int i = 0; i < size; i++ )
          {
            const double bl = 0.3;
            double vn = ( 1.0 - bl ) * vi[0][i] + bl * vo[0][i];
            vi1[0][i] = vi[0][i];
            vo1[0][i] = vo[0][i];
            vi[0][i]  = vn;
            
          }
        }
        else
        {
          // itscf > 1
          double sn = 0.0;
          double sd = 0.0;
          for ( int i = 0; i < size; i++ )
          {
            double rl = vo[0][i] - vi[0][i];
            double rl1 = vo1[0][i] - vi1[0][i];
            double dr = rl - rl1;
            sn += rl * dr;
            sd += dr * dr;
          }
          double tmp[2];
          tmp[0] = sn;
          tmp[1] = sd;
          s_.wf.context().dsum(2,1,&tmp[0],2);
          sn = tmp[0];
          sd = tmp[1];
          
          if ( sd != 0.0 )
            thl = sn / sd;
          else
            thl = 1.0;

          if ( s_.ctxt_.onpe0() )
            cout << "  <!-- potential mixing: Anderson thl=" << thl << " -->"
                 << endl;
          for ( int i = 0; i < size; i++ )
          {
            const double bl = 0.3;
            double vn = ( 1.0 - bl ) * ( ( 1.0 - thl ) * 
                        vi[0][i] + thl * vi1[0][i] ) +
                        bl * ( ( 1.0 - thl ) * vo[0][i] + thl * vo1[0][i] );
            vi1[0][i] = vi[0][i];
            vo1[0][i] = vo[0][i];
            vi[0][i]  = vn;
          }
        }
        ef_.v_r = vi;
        
        }
#endif

        // Next line: reset the wf stepper only if nite > 1
        // If nite = 1, the acceleration procedure of some steppers should not 
        // be reset.
        if ( nite_ > 1 ) wf_stepper->preprocess();
      
        for ( int ite = 0; ite < nite_; ite++ )
        {
          double energy = 0.0;

          if ( compute_forces && compute_fv )
          {
            // compute forces at each electronic step for monitoring
            energy = ef_.energy(true,dwf,compute_forces,fion,false,sigma_eks);
          }
          else
          {
            // normal case: do not compute forces during wf optimization
            energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
          }
          
          // compute the sum of eigenvalues (with fixed weight)
          // to measure convergence of the subspace update
          // compute trace of the Hamiltonian matrix Y^T H Y
          // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
          const double eigenvalue_sum = s_.wf.dot(dwf);
          if ( s_.ctxt_.onpe0() )
            cout << "  <eigenvalue_sum> "
                 << eigenvalue_sum << " </eigenvalue_sum>" << endl;
 
          wf_stepper->update(dwf);
          
          if ( s_.ctxt_.onpe0() )
          {
            cout.setf(ios::fixed,ios::floatfield);
            cout.setf(ios::right,ios::adjustfield);
//             cout << "  <ekin_int>   " << setprecision(8)
//                  << setw(15) << ef_.ekin() << " </ekin_int>\n";
//             if ( use_confinement )
//             {
//               cout << "  <econf_int>  " << setw(15) << ef_.econf()
//                    << " </econf_int>\n";
//             }
//             cout << "  <eps_int>    " << setw(15) 
//                  << ef_.eps() << " </eps_int>\n"
//                  << "  <enl_int>    " << setw(15) 
//                  << ef_.enl() << " </enl_int>\n"
//                  << "  <ecoul_int>  " << setw(15) 
//                  << ef_.ecoul() << " </ecoul_int>\n"
//                  << "  <exc_int>    " << setw(15) 
//                  << ef_.exc() << " </exc_int>\n"
//                  << "  <esr_int>    " << setw(15) 
//                  << ef_.esr() << " </esr_int>\n"
//                  << "  <eself_int>  " << setw(15) 
//                  << ef_.eself() << " </eself_int>\n"
              cout << "  <etotal_int> " << setw(15) 
                   << ef_.etotal() << " </etotal_int>\n";
            if ( compute_stress )
            {
              const double pext = (sigma_ext[0]+sigma_ext[1]+sigma_ext[2])/3.0;
              const double enthalpy = ef_.etotal() + pext * cell.volume();
              cout << "  <enthalpy_int> " << setw(15) 
                   << enthalpy << " </enthalpy_int>\n"
                   << flush;
            }
 
            // compute force*velocity
            if ( compute_forces && compute_fv )
            {
              double fv = 0.0;
              for ( int is = 0; is < atoms.atom_list.size(); is++ )
              {
                int i = 0;
                for ( int ia = 0; ia < atoms.atom_list[is].size(); ia++ )
                {
                  fv += ionic_stepper->v0(is,i) * fion[is][i] +
                        ionic_stepper->v0(is,i+1) * fion[is][i+1] +
                        ionic_stepper->v0(is,i+2) * fion[is][i+2];
                  i += 3;
                }
              }
              cout << "  <fv> " << fv << " </fv>\n";
            }
          }
        } // for ite
        
        // subspace diagonalization
        if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" )
        {
          energy = ef_.energy(true,dwf,false,fion,false,sigma_eks);
          s_.wf.diag(dwf,compute_eigvec);
        }
        
        // update occupation numbers
        if ( nempty > 0 )
        {
          s_.wf.update_occ(s_.ctrl.fermi_temp);
          const double wf_entropy = s_.wf.entropy();
          if ( s_.ctxt_.onpe0() )
          {
            cout << "  <!-- Wavefunction entropy: " << wf_entropy
                 << " -->" << endl;
            const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
            cout << "  <!-- Entropy contribution to free energy: "
                 << - wf_entropy * s_.ctrl.fermi_temp * boltz
                 << " -->" << endl;
          }
        }
        
        if ( s_.ctxt_.onpe0() )
          cout << "  <!-- BOSampleStepper: end scf iteration -->" << endl;
      } // for itscf
      
      wf_stepper->postprocess();
    }
    else
    {
      // wf_stepper == 0, wf_dyn == LOCKED
      // evaluate and print energy
      cd_.update_density();
      ef_.update_vhxc();
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
    cd_.update_density();
    ef_.update_vhxc();
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
