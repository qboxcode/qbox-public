////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.C,v 1.16 2003-10-02 17:41:04 fgygi Exp $

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleStepper::SampleStepper(Sample& s) : s_(s), dwf(s.wf), wfv(s.wfv)
{
  fion.resize(s_.atoms.nsp());
  tau0.resize(s_.atoms.nsp());
  taum.resize(s_.atoms.nsp());
  vel.resize(s_.atoms.nsp());
  pmass.resize(s_.atoms.nsp());
  for ( int is = 0; is < fion.size(); is++ )
  {
    fion[is].resize(3*s_.atoms.na(is));
    tau0[is].resize(3*s_.atoms.na(is));
    taum[is].resize(3*s_.atoms.na(is));
    vel[is].resize(3*s_.atoms.na(is));
    pmass[is] = s_.atoms.species_list[is]->mass() * 1822.89;
  }
}

////////////////////////////////////////////////////////////////////////////////
SampleStepper::~SampleStepper(void)
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
void SampleStepper::step(EnergyFunctional& e, int niter)
{
  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  atoms.get_positions(tau0);
  atoms.get_velocities(vel);
  
  const double dt = s_.ctrl.dt;
  // default value of emass if 
  const double emass = s_.ctrl.emass;
  const double dt2bye = (emass == 0.0) ? 0.5 / wf.ecut() : dt*dt/emass;            
  double ekin_ion, ekin_e, ekin_em=0.0, temp_ion, eta;
  
  const bool compute_hpsi = ( s_.ctrl.wf_dyn != "LOCKED" );
  const bool compute_forces = ( s_.ctrl.atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool gram_ortho = ( s_.ctrl.wf_dyn != "MD" );
  const bool ttherm = ( s_.ctrl.thermostat == "ON" );
  const int ndofs = 3 * s_.atoms.size();
  const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
  const double th_temp = s_.ctrl.th_temp;
  const double th_time = s_.ctrl.th_time;
  
  const bool onrow0 = ( wf.context().myrow() == 0 );
  
  const double ecutprec = s_.ctrl.ecutprec;
  const bool precondition = ( s_.ctrl.wf_dyn == "PSD" && ecutprec != 0.0 );
  
  if ( s_.ctxt_.myproc()==0 )
  {
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        if ( wf.sd(ispin,ikp) != 0 )
        {
          if ( wf.sdcontext(ispin,ikp)->active() )
          {
            SlaterDet* sd = wf.sd(ispin,ikp);
            cout << "  <!-- kpoint = " << sd->kpoint() << " -->" << endl;
            cout << "  <!-- basis.size() = " << sd->basis().size() 
                 << " -->" << endl;
            cout << "  <!-- nst = " << sd->nst() << " -->" << endl;
            cout << "  <!-- c dimensions: "
                 << sd->c().m() << "x" << sd->c().n() 
                 << "   (" << sd->c().mb() << "x"
                 << sd->c().nb() << " blocks)" << " -->" << endl;
          }
        }
      }
    }
  }
  
  Timer tm_iter;
 
  for ( int iter = 0; iter < niter; iter++ )
  {
  
  tm_iter.start();
  
  if ( s_.ctxt_.mype() == 0 )
    cout << "  <iteration count=\"" << iter+1 << "\">\n";
    
  double energy =
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
         << "  <etotal> " << setw(15) << e.etotal() << " </etotal>\n" << flush;
  }
  
  if ( compute_hpsi )
  {
    double ekin_ep = 0.0;
    double ekin_e0 = 0.0;
    // update SlaterDet
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        if ( wf.sd(ispin,ikp) != 0 )
        {
          if ( wf.sdcontext(ispin,ikp)->active() )
          {
            if ( s_.ctrl.wf_diag == "T" )
            {
            tmap["eigval"].start();
            // compute eigenvalues
            if ( wf.sd(ispin,ikp)->basis().real() )
            {
              // proxy real matrices c, cp
              DoubleMatrix c(wf.sd(ispin,ikp)->c());
              DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
              DoubleMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
              // factor 2.0 in next line: G and -G
              h.gemm('t','n',2.0,c,cp,0.0);
              // rank-1 update correction
              h.ger(-1.0,c,0,cp,0);
 
              //cout << " Hamiltonian at k = " << wf.sd(ispin,ikp)->kpoint()
              //     << endl;
              //cout << h;
 
              valarray<double> w(h.m());
              h.syev('l',w);
              if ( s_.wf.context().onpe0() )
              {
                const double eVolt = 2.0 * 13.6058;
                cout <<    "  <eigenvalues spin=\"" << ispin
                     << "\" kpoint=\"" << wf.sd(ispin,ikp)->kpoint()
                     << "\" n=\"" << h.m() << "\">" << endl;
                for ( int i = 0; i < h.m(); i++ )
                {
                  cout << setw(10) << setprecision(5) << w[i]*eVolt;
                  if ( i%5 == 4 ) cout << endl;
                }
                if ( h.m()%5 != 0 ) cout << endl;
                cout << "  </eigenvalues>" << endl;
              }
            }
            else
            {
              // complex case not implemented
              assert(false);
              #if 0
              ComplexMatrix& c(wf.sd[ikp]->c());
              ComplexMatrix& cp(dwf.sd[ikp]->c());
 
              ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
              h.gemm('c','n',1.0,c,cp,0.0);
 
              //cout << " Hamiltonian at k = " << wf.sd[ikp]->kpoint() << endl;
              //cout << h;
 
              valarray<double> w(h.m());

              h.heev('l',w);
              cout << " Eigenvalues at k = " << wf.sd[ikp]->kpoint() << endl;
              const double eVolt = 2.0 * 13.6058;
              for ( int i = 0; i < h.m(); i++ )
              {
                cout << "%" << setw(3) << ikp
                     << setw(10) << setprecision(5) << w[i]*eVolt << endl;;
              }
              #endif
            }
            tmap["eigval"].stop();
 
            } // wfdiag T

            tmap["update_psi"].start();
            
            if ( precondition )
            {
              // compute A = V^T H V  and descent direction HV - VA
              
              if ( wf.sd(ispin,ikp)->basis().real() )
              {
                // proxy real matrices c, cp
                DoubleMatrix c(wf.sd(ispin,ikp)->c());
                DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
                DoubleMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
                // factor 2.0 in next line: G and -G
                a.gemm('t','n',2.0,c,cp,0.0);
                // rank-1 update correction
                a.ger(-1.0,c,0,cp,0);
                
                // cp = cp - c * a
                cp.gemm('n','n',-1.0,c,a,1.0);
              }
              else
              {
                // not implemented in the complex case
                assert(false);
              }
              
              // dwf.sd->c() now contains the descent direction (HV-VA)
 
              const double g2i_prec = 0.5 / ecutprec;
              const double* g2i_ptr = wf.sd(ispin,ikp)->basis().g2i_ptr();
              double* coeff = (double*) wf.sd(ispin,ikp)->c().valptr();
              const double* dcoeff = 
                (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
              const int mloc = wf.sd(ispin,ikp)->c().mloc();
              const int nloc = wf.sd(ispin,ikp)->c().nloc();
              for ( int n = 0; n < nloc; n++ )
              {
                // note: double mloc length for complex<double> indices
                double* c = &coeff[2*mloc*n];
                const double* dc = &dcoeff[2*mloc*n];
                for ( int i = 0; i < mloc; i++ )
                {
                  const double g2i = g2i_ptr[i];
                  const double dt2bye = ( g2i == 0.0 ? g2i_prec : 
                    ( g2i < g2i_prec ) ? g2i : g2i_prec );
                  c[2*i] -= dt2bye * dc[2*i];
                  c[2*i+1] -= dt2bye * dc[2*i+1];
                }
              }
            }
            else
            {
              // no preconditioning
              // update wf using SD or MD
              
              if ( s_.ctrl.wf_dyn == "MD" )
              {
                assert(wfv!=0);
                if ( iter == 0 )
                {
                  // First iteration of Stoermer's rule
                  // compute wf from velocity wfv
                  // cp = c + dt * v - 0.5 * dt2/m * hpsi
                  // cm = c
                  // ekin_e0 += 0.5 * m_e * v * v
                  // 0.5 * m_e = 0.5 * dt * dt / dt2bye
                  const double half_m_e = 0.5 * dt * dt / dt2bye;
                  const double half_dt2bye = 0.5 * dt2bye;
                  double* cptr = (double*) wf.sd(ispin,ikp)->c().valptr();
                  double* cptrv = (double*) wfv->sd(ispin,ikp)->c().valptr();
                  const double* dcptr =
                    (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
                  const vector<double>& occ = wf.occ(ispin,ikp);
                  const int mloc = wf.sd(ispin,ikp)->c().mloc();
                  const int nloc = wf.sd(ispin,ikp)->c().nloc();
                  for ( int n = 0; n < nloc; n++ )
                  {
                    const int nglobal = wf.sd(ispin,ikp)->c().j(0,n);
                    const double occn = occ[nglobal];
                    // note: double mloc length for complex<double> indices
                    double* c = &cptr[2*mloc*n];
                    double* cv = &cptrv[2*mloc*n];
                    const double* dc = &dcptr[2*mloc*n];
                    double tmpsum = 0.0;
                    if ( onrow0 )
                    {
                      const double cvtmp = cv[0];
                      tmpsum -= 0.5 * (cvtmp*cvtmp);
                      cv[1] = 0.0;
                    }
                    for ( int i = 0; i < mloc; i++ )
                    {
                      const double ctmp = c[2*i];
                      const double ctmp1 = c[2*i+1];
                      const double cvtmp = cv[2*i];
                      const double cvtmp1 = cv[2*i+1];
                      const double dctmp = dc[2*i];
                      const double dctmp1 = dc[2*i+1];
 
                      tmpsum += (cvtmp*cvtmp + cvtmp1*cvtmp1);
                      c[2*i]    = ctmp  + dt * cvtmp  - half_dt2bye * dctmp;
                      c[2*i+1]  = ctmp1 + dt * cvtmp1 - half_dt2bye * dctmp1;
                      cv[2*i]   = ctmp;
                      cv[2*i+1] = ctmp1;
                    }
                    // Note: 2 in next line: from G,-G
                    //!! extra factor of 2
                    ekin_e0 += 4.0 * half_m_e * occn * tmpsum;
                  }
                  // Note: *wfv now contains wf(t-dt)
                }
                else
                {
                  // Verlet update of wf
                  // cp = c + (c - cm) + dt2/m * hpsi
                  // c += c - cm + dt2bye * hpsi
                  // cm = c
                  double* cptr = (double*) wf.sd(ispin,ikp)->c().valptr();
                  double* cptrm = (double*) wfv->sd(ispin,ikp)->c().valptr();
                  const double* dcptr =
                    (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
                  const int mloc = wf.sd(ispin,ikp)->c().mloc();
                  const int nloc = wf.sd(ispin,ikp)->c().nloc();
                  for ( int n = 0; n < nloc; n++ )
                  {
                    // note: double mloc length for complex<double> indices
                    double* c = &cptr[2*mloc*n];
                    double* cm = &cptrm[2*mloc*n];
                    const double* dc = &dcptr[2*mloc*n];
                    for ( int i = 0; i < mloc; i++ )
                    {
                      const double ctmp = c[2*i];
                      const double ctmp1 = c[2*i+1];
                      const double cmtmp = cm[2*i];
                      const double cmtmp1 = cm[2*i+1];
                      const double dctmp = dc[2*i];
                      const double dctmp1 = dc[2*i+1];
                      const double cptmp = 2.0*ctmp -  cmtmp -  dt2bye * dctmp;
                      const double cptmp1= 2.0*ctmp1 - cmtmp1 - dt2bye * dctmp1;
 
                      c[2*i]    = cptmp;
                      c[2*i+1]  = cptmp1;
                      cm[2*i]   = ctmp;
                      cm[2*i+1] = ctmp1;
                    }
                  }
                }
              }
              else
              {
                // SD update of wf
                // c = c - dt2bye * hpsi
                wf.sd(ispin,ikp)->c().axpy(-dt2bye,dwf.sd(ispin,ikp)->c());
              }
            }
            tmap["update_psi"].stop();

            if ( gram_ortho )
            {
              tmap["gram"].start();
              wf.sd(ispin,ikp)->gram();
              tmap["gram"].stop();
            }
            else
            {
              tmap["riccati"].start();
              assert(wfv!=0);
              wf.sd(ispin,ikp)->riccati(*wfv->sd(ispin,ikp));
              tmap["riccati"].stop();
              
              // compute electronic kinetic energy at time t+1/2
              double* cptr = (double*) wf.sd(ispin,ikp)->c().valptr();
              double* cptrm = (double*) wfv->sd(ispin,ikp)->c().valptr();
              const vector<double>& occ = wf.occ(ispin,ikp);
              const int mloc = wf.sd(ispin,ikp)->c().mloc();
              const int nloc = wf.sd(ispin,ikp)->c().nloc();
              for ( int n = 0; n < nloc; n++ )
              {
                const int nglobal = wf.sd(ispin,ikp)->c().j(0,n);
                const double occn = occ[nglobal];
                // note: double mloc length for complex<double> indices
                double* c = &cptr[2*mloc*n];
                double* cm = &cptrm[2*mloc*n];
                double tmpsum = 0.0;
                if ( onrow0 )
                {
                  // correct for double counting of G=0 element
                  // i=0 coefficient is real, count only real part
                  const double ctmp = c[0];
                  const double cmtmp = cm[0];
                  tmpsum -= 0.5 * (ctmp - cmtmp)*(ctmp - cmtmp);
                }
                for ( int i = 0; i < mloc; i++ )
                {
                  const double ctmp = c[2*i];
                  const double ctmp1 = c[2*i+1];
                  const double cmtmp = cm[2*i];
                  const double cmtmp1 = cm[2*i+1];
 
                  tmpsum += (ctmp -cmtmp )*(ctmp -cmtmp ) +
                            (ctmp1-cmtmp1)*(ctmp1-cmtmp1);
                }
                // Note: 2 in next line: from (G,-G)
                //!! factor 1/2 from 1/2 m v^2 absent
                ekin_ep += ( 2.0 * occn / dt2bye ) * tmpsum;
              }
            }
          }
        }
      }
    }
    if ( iter == 0 )
    {
      ekin_e = ekin_e0;
    }
    else
    {
      ekin_e = 0.5 * ( ekin_em + ekin_ep );
    }
    ekin_em = ekin_ep;
    wf.context().dsum(1,1,&ekin_e,1);
    if ( s_.wf.context().onpe0() && s_.ctrl.wf_dyn == "MD" )
    {
      cout << "  <ekin_e> " << ekin_e << " </ekin_e>\n";
    }
  } // if compute_hpsi
  
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
               << "    <position> " << tau0[is][i] << " "
               << tau0[is][i+1] << " " << tau0[is][i+2] << " </position>\n"
               << "    <force> " << fion[is][i] << " "
               << fion[is][i+1] << " " << fion[is][i+2] 
               << " </force>\n  </atom>" << endl;
               
          i += 3;
        }
      }
    }
    
    if ( s_.ctrl.atoms_dyn == "SD" )
    {
      // update tau0 (SD)
      for ( int is = 0; is < tau0.size(); is++ )
      {
        const double dt2bym = dt * dt / pmass[is];
        for ( int i = 0; i < tau0[is].size(); i++ )
        {
          tau0[is][i] += dt2bym * fion[is][i];
        }
      }
    }
    else
    {
      // MD for ions
      if ( iter == 0 )
      {
        // First step of Stoermer's rule
        // x1 = x0 + h * v0 + 0.5 * dt2/m * f0
        
        // compute ekin_ion and temp_ion
        ekin_ion = 0.0;
        eta = 0.0;
        for ( int is = 0; is < tau0.size(); is++ )
        {
          const double dt2bym = dt * dt / pmass[is];
          for ( int i = 0; i < tau0[is].size(); i++ )
          {
            ekin_ion += 0.5 * pmass[is] * vel[is][i] * vel[is][i];
          }
        }
        
        if ( ndofs > 0 )
          temp_ion = 2.0 * ( ekin_ion / boltz ) / ndofs;
        else
          temp_ion = 0.0;

        if ( ttherm )
        {
          // Next line: linear thermostat, width of tanh = 100.0
          eta = tanh ( ( temp_ion - th_temp ) / 100. ) / th_time;
        }
        
        for ( int is = 0; is < tau0.size(); is++ )
        {
          const double dt2bym = dt * dt / pmass[is];
          for ( int i = 0; i < tau0[is].size(); i++ )
          {
            const double taup = tau0[is][i] + dt * vel[is][i] -
            dt * dt * eta * vel[is][i] + 
            0.5 * dt2bym * fion[is][i];
            taum[is][i] = tau0[is][i];
            tau0[is][i] = taup;
          }
        }
      }
      else
      {
        // Normal Verlet step
        // compute ekin_ion and temp_ion before step using a first order
        // approximation. Note: ekin_ion is recomputed after the step using
        // a second-order approximation.
        ekin_ion = 0.0;
        eta = 0.0;
        for ( int is = 0; is < tau0.size(); is++ )
        {
          const double dt2bym = dt * dt / pmass[is];
          if ( dt != 0.0 )
          {
            for ( int i = 0; i < tau0[is].size(); i++ )
            {
              const double v = ( tau0[is][i] - taum[is][i] ) / dt;
              ekin_ion += 0.5 * pmass[is] * v * v;
            }
          }
        }
        
        if ( ndofs > 0 )
          temp_ion = 2.0 * ( ekin_ion / boltz ) / ndofs;
        else
          temp_ion = 0.0;

        if ( ttherm )
        {
          // Next line: linear thermostat, width of tanh = 100.0
          eta = tanh ( ( temp_ion - th_temp ) / 100. ) / th_time;
        }
        
        ekin_ion = 0.0;
        for ( int is = 0; is < tau0.size(); is++ )
        {
          const double dt2bym = dt *dt / pmass[is];
          for ( int i = 0; i < tau0[is].size(); i++ )
          {
            const double taup = tau0[is][i] +
            (1.0 - dt*eta) * (tau0[is][i] - taum[is][i]) +
            dt2bym * fion[is][i];
            if ( dt != 0.0 )
            {
              const double v = 0.5 * ( taup - taum[is][i] ) / dt;
              ekin_ion += 0.5 * pmass[is] * v * v;
            }
            taum[is][i] = tau0[is][i];
            tau0[is][i] = taup;
          }
        }
        
        if ( ndofs > 0 )
          temp_ion = 2.0 * ( ekin_ion / boltz ) / ndofs;
        else
          temp_ion = 0.0;

      }
      if ( s_.wf.context().onpe0() )
      {
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
        cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
        cout << "  <eta_ion> " << eta << " </eta_ion>\n";
      }
    }
    atoms.set_positions(tau0);
    e.atoms_moved();
  }
  if ( s_.wf.context().onpe0() )
  {
    cout << "  <econst> " << energy+ekin_ion+ekin_e << " </econst>\n";
    if ( s_.ctrl.wf_dyn == "MD" )
      cout << "  <ekin_ec> " << energy+ekin_ion+2*ekin_e << " </ekin_ec>\n";
  }
 
  if ( compute_stress )
  {
    // update unit cell
    e.cell_moved();
  }
  
  // print iteration time
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
    
  if ( s_.ctxt_.mype() == 0 )
    cout << "  </iteration>" << endl;
  } // niter
  
  // if MD was used, compute final velocity of ions and wavefunctions
  // The calculation of the final velocity require forces at the endpoint
  if ( s_.ctrl.wf_dyn == "MD" || s_.ctrl.atoms_dyn == "MD" )
  {
    double energy =
       e.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,dcell);
  }
  // dwf and fion now contain the forces on wavefunctions and ions at the 
  // endpoint
  
  // final velocity of ions
  if ( compute_forces )
  {
    if ( s_.ctrl.atoms_dyn == "MD" )
    {
      // Last step of Stoermer's rule
      // v = (xn - x(n-1) )/dt + 0.5 * dt/m * fn
         
      if ( dt != 0.0 )
      {
        for ( int is = 0; is < tau0.size(); is++ )
        {
          const double dtbym = dt / pmass[is];
          for ( int i = 0; i < tau0[is].size(); i++ )
          {
            vel[is][i] = ( tau0[is][i] - taum[is][i] ) / dt +
                         0.5 * dtbym * fion[is][i];
          }
        }
      }
    }
    else
    {
      // for all methods other than MD, set velocity to zero
      for ( int is = 0; is < tau0.size(); is++ )
      {
        const double dtbym = dt / pmass[is];
        for ( int i = 0; i < tau0[is].size(); i++ )
        {
          vel[is][i] = 0.0;
        }
      }
    }
    atoms.set_velocities(vel);
  }
  
  // final velocity of wavefunctions
  if ( compute_hpsi )
  {
    if ( s_.ctrl.wf_dyn == "MD" )
    {
      assert(wfv!=0);
      if ( dt != 0.0 )
      {
        for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
        {
          for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
          {
            if ( wf.sd(ispin,ikp) != 0 )
            {
              if ( wf.sdcontext(ispin,ikp)->active() )
              {
                // Last iteration of Stoermer's rule
                // compute final velocity wfv
                // v = ( c - cm ) / dt - 0.5 * dt/m * hpsi
 
                // Note: At this point, *wfv contains wf(t-dt)
                
                // hpsi must be orthogonal to subspace spanned by c
                // compute descent direction H psi - psi (psi^T H psi)
                
                if ( wf.sd(ispin,ikp)->basis().real() )
                {
                  // proxy real matrices c, cp
                  DoubleMatrix c(wf.sd(ispin,ikp)->c());
                  DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
                  DoubleMatrix a(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
                  // factor 2.0 in next line: G and -G
                  a.gemm('t','n',2.0,c,cp,0.0);
                  // rank-1 update correction
                  a.ger(-1.0,c,0,cp,0);
 
                  // cp = cp - c * a
                  cp.gemm('n','n',-1.0,c,a,1.0);
                }
                else
                {
                  // not implemented in the complex case
                  assert(false);
                }
              
                const double dt_inv = 1.0/dt;
                const double half_dtbye = 0.5 * dt2bye / dt;
                double* cptr = (double*) wf.sd(ispin,ikp)->c().valptr();
                double* cptrv = (double*) wfv->sd(ispin,ikp)->c().valptr();
                const double* dcptr =
                  (const double*) dwf.sd(ispin,ikp)->c().cvalptr();
                const int mloc = wf.sd(ispin,ikp)->c().mloc();
                const int nloc = wf.sd(ispin,ikp)->c().nloc();
                for ( int n = 0; n < nloc; n++ )
                {
                  // note: double mloc length for complex<double> indices
                  double* c = &cptr[2*mloc*n];
                  double* cv = &cptrv[2*mloc*n];
                  const double* dc = &dcptr[2*mloc*n];
                  for ( int i = 0; i < mloc; i++ )
                  {
                    const double ctmp = c[2*i];
                    const double ctmp1 = c[2*i+1];
                    const double cmtmp = cv[2*i];
                    const double cmtmp1 = cv[2*i+1];
                    const double dctmp = dc[2*i];
                    const double dctmp1 = dc[2*i+1];
 
                    cv[2*i]   = ( ctmp  - cmtmp  ) * dt_inv 
                                - half_dtbye * dctmp;
                    cv[2*i+1] = ( ctmp1 - cmtmp1 ) * dt_inv 
                                - half_dtbye * dctmp1;
                  }
                }
                // Note: *wfv now contains the wavefunction velocity
              }
            }
          }
        }
      }
    }
  }
}
