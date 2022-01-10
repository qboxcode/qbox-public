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
// BOSampleStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "BOSampleStepper.h"
#include "EnergyFunctional.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "WavefunctionStepper.h"
#include "SDWavefunctionStepper.h"
#include "PSDWavefunctionStepper.h"
#include "PSDAWavefunctionStepper.h"
#include "JDWavefunctionStepper.h"
#include "UserInterface.h"
#include "Preconditioner.h"
#include "SDIonicStepper.h"
#include "SDAIonicStepper.h"
#include "CGIonicStepper.h"
#include "ANDIonicStepper.h"
#include "MDIonicStepper.h"
#include "BMDIonicStepper.h"
#include "SDCellStepper.h"
#include "CGCellStepper.h"
#include "AndersonMixer.h"
#include "MLWFTransform.h"
#include "D3tensor.h"
#include "cout0.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::BOSampleStepper(Sample& s, int nitscf, int nite) :
  SampleStepper(s), cd_(s.wf), ef_(s,cd_),
  dwf(s.wf), nitscf_(nitscf), nite_(nite),
  update_density_first_(true), update_vh_(true), update_vxc_(true) {}

////////////////////////////////////////////////////////////////////////////////
BOSampleStepper::~BOSampleStepper()
{
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
void BOSampleStepper::initialize_density(void)
{
  // initialize cd_ with a sum of atomic densities

  double atom_radius[] =
  {
    0.000, 1.542, 0.931, 1.727, 1.636, // null H  He  Li  Be
    1.845, 1.538, 1.311, 1.141, 1.014, //  B   C   N   O   F
    0.983, 1.238, 1.347, 1.376, 2.151, // Ne  Na  Mg  Al  Si
    1.927, 1.733, 1.563, 1.429, 1.546, //  P   S  Cl  Ar   K
    1.663, 1.604, 1.516, 1.434, 1.377, // Ca  Sc  Ti   V  Cr
    1.310, 1.245, 1.201, 1.162, 1.098, // Mn  Fe  Co  Ni  Cu
    1.077, 1.331, 1.415, 2.015, 1.880, // Zn  Ga  Ge  As  Se
    1.749, 1.630, 1.705, 1.819, 1.794, // Br  Kr  Rb  Sr   Y
    1.728, 1.664, 1.589, 1.523, 1.461, // Zr  Nb  Mo  Tc  Ru
    1.410, 1.348, 1.306, 1.303, 1.554, // Rh  Pd  Ag  Cd  In
    1.609, 1.611, 1.530, 1.514, 1.464, // Sn  Sb  Te   I  Xe
    1.946, 1.967, 1.943, 1.930, 1.920, // Cs  Ba  La  Ce  Pr
    1.910, 1.900, 1.890, 1.880, 1.870, // Nd  Pm  Sm  Eu  Gd
    1.860, 1.850, 1.840, 1.830, 1.820, // Tb  Dy  Ho  Er  Tm
    1.810, 1.800, 1.701, 1.658, 1.606, // Yb  Lu  Hf  Ta   W
    1.550, 1.500, 1.446, 1.398, 1.355, // Re  Os  Ir  Pt  Au
    1.314, 1.624, 1.659, 1.634, 1.620, // Hg  Tl  Pb  Bi  Po
    1.600, 1.500, 1.600, 1.600, 1.600, // At  Rn  Fr  Ra  Ac
    1.600, 1.600, 1.600, 1.600, 1.600, // Th  Pa   U  Np  Pu
    1.600, 1.600, 1.600, 1.600, 1.600, // Am  Cm  Bk  Cf  Es
    1.600, 1.600, 1.600, 1.600         // Fm  Md  No  Lr
  };

  const AtomSet& atoms = s_.atoms;
  const Basis* const vbasis = cd_.vbasis();
  const int ngloc = vbasis->localsize();
  vector<vector<complex<double> > > rhops;
  const int nsp = atoms.nsp();
  rhops.resize(nsp);
  vector<complex<double> > rhopst(ngloc);
  const double * const g2 = vbasis->g2_ptr();

  for ( int is = 0; is < nsp; is++ )
  {
    rhops[is].resize(ngloc);
    Species *s = atoms.species_list[is];
    const int zval = s->zval();
    const int atomic_number = s->atomic_number();
    assert(atomic_number < sizeof(atom_radius)/sizeof(double));
    double rc = atom_radius[atomic_number];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      double arg = 0.25 * rc * rc * g2[ig];
      rhops[is][ig] = zval * exp( -arg );
    }
  }

  vector<vector<double> > tau0;
  tau0.resize(nsp);
  for ( int is = 0; is < nsp; is++ )
    tau0.resize(3*atoms.na(is));
  atoms.get_positions(tau0);
  StructureFactor sf;
  sf.init(tau0,*vbasis);
  sf.update(tau0,*vbasis);

  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < nsp; is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
    }
  }

  // Adjust G=0 component of the charge if net_charge is non-zero
  rhopst[0] += s_.wf.nel() - atoms.nel();

  // Initialize charge equally for both spins
  cd_.rhog[0] = rhopst;
  if ( cd_.rhog.size() == 2 )
  {
    assert(cd_.rhog[0].size()==cd_.rhog[1].size());
    for ( int i = 0; i < cd_.rhog[0].size(); i++ )
    {
      cd_.rhog[0][i] = 0.5 * rhopst[i];
      cd_.rhog[1][i] = 0.5 * rhopst[i];
    }
  }
  cd_.update_rhor();
  update_density_first_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void BOSampleStepper::step(int niter)
{
  const bool anderson_charge_mixing = ( s_.ctrl.charge_mix_ndim > 0 );

  // determine whether eigenvectors must be computed
  // eigenvectors are computed if explicitly requested with wf_diag==T
  // or if the SlaterDet has fractionally occupied states
  const bool fractional_occ = (s_.wf.nel() != 2 * s_.wf.nst());
  const bool compute_eigvec = fractional_occ || s_.ctrl.wf_diag == "T";
  const bool compute_mlwf = s_.ctrl.wf_diag == "MLWF";
  const bool compute_mlwfc = s_.ctrl.wf_diag == "MLWFC";
  enum ortho_type { GRAM, LOWDIN, ORTHO_ALIGN, RICCATI };

  const double gpa = 29421.5;

  AtomSet& atoms = s_.atoms;
  Wavefunction& wf = s_.wf;
  Wavefunction*& wfv = s_.wfv;
  const bool onpe0 = MPIdata::onpe0();
  const Context& sd_ctxt = wf.sd_context();

  const int nspin = wf.nspin();

  const double dt = s_.ctrl.dt;

  const string wf_dyn = s_.ctrl.wf_dyn;
  const string atoms_dyn = s_.ctrl.atoms_dyn;
  const string cell_dyn = s_.ctrl.cell_dyn;

  const bool extrapolate_wf = ( atoms_dyn == "MD" );

  bool ntc_extrapolation = false;
  bool asp_extrapolation = false;

  const map<string,string>& debug_map = s_.ctrl.debug;

  map<string,string>::const_iterator imap =
    debug_map.find("EXTRAPOLATION");
  if ( imap != debug_map.end() )
  {
    const string val = imap->second;
    if ( val == "NTC" ) ntc_extrapolation = true;
    if ( val == "ASP" ) asp_extrapolation = true;
  }

  Wavefunction* wfmm;
  if ( extrapolate_wf && ( ntc_extrapolation || asp_extrapolation ) )
    wfmm = new Wavefunction(wf);

  // Next lines: special value of niter = 0: GS calculation only
  const bool atoms_move = ( niter > 0 && atoms_dyn != "LOCKED" );
  const bool compute_stress = ( s_.ctrl.stress == "ON" );
  const bool cell_moves = ( niter > 0 && compute_stress &&
                            cell_dyn != "LOCKED" );
  // GS-only calculation:
  const bool gs_only = !atoms_move && !cell_moves;

  const double force_tol = s_.ctrl.force_tol;
  const double stress_tol = s_.ctrl.stress_tol;

  Timer tm_iter;

  Preconditioner prec(wf,ef_,s_.ctrl.ecutprec);

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
    wf_stepper = new PSDWavefunctionStepper(wf,prec,tmap);
  else if ( wf_dyn == "PSDA" )
    wf_stepper = new PSDAWavefunctionStepper(wf,prec,tmap);
  else if ( wf_dyn == "JD" )
    wf_stepper = new JDWavefunctionStepper(wf,prec,ef_,tmap);

  // wf_stepper == 0 indicates that wf_dyn == LOCKED

  IonicStepper* ionic_stepper = 0;
  if ( atoms_dyn == "SD" )
    ionic_stepper = new SDIonicStepper(s_);
  else if ( atoms_dyn == "SDA" )
    ionic_stepper = new SDAIonicStepper(s_);
  else if ( atoms_dyn == "CG" )
    ionic_stepper = new CGIonicStepper(s_);
  else if ( atoms_dyn == "AND" )
    ionic_stepper = new ANDIonicStepper(s_);
  else if ( atoms_dyn == "MD" )
    ionic_stepper = new MDIonicStepper(s_);
  else if ( atoms_dyn == "BMD" )
    ionic_stepper = new BMDIonicStepper(s_);

  if ( ionic_stepper )
    ionic_stepper->setup_constraints();

  CellStepper* cell_stepper = 0;
  if ( cell_dyn == "SD" )
    cell_stepper = new SDCellStepper(s_);
  else if ( cell_dyn == "CG" )
    cell_stepper = new CGCellStepper(s_);

  // Allocate wavefunction velocity if not available
  if ( atoms_move && extrapolate_wf )
  {
    if ( wfv == 0 )
    {
      wfv = new Wavefunction(wf);
      wfv->clear();
    }
  }

  vector<MLWFTransform*> mlwft(wf.nsp_loc());

  if ( compute_mlwf || compute_mlwfc )
  {
    // MLWF can be computed at the gamma point only
    // There must be a single k-point, and it must be gamma
    if ( wf.nkp() > 1 || ( wf.nkp()==1 && wf.kpoint(0) != D3vector(0,0,0) ) )
    {
      if ( onpe0 )
      {
        cout << " BOSampleStepper::step: MLWF can be computed at k=0 only"
             << endl;
        cout << " BOSampleStepper::step: cannot run" << endl;
      }
      return;
    }

    for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
    {
      const int ikp = 0;
      const int ikp_loc = wf.ikp_local(ikp);
      if ( ikp_loc >= 0 )
        mlwft[isp_loc] = new MLWFTransform(*wf.sd(isp_loc,ikp_loc));
    }
  }

  // Charge mixing variables: include both spins in the same vector
  if ( nspin > 1 ) assert(cd_.rhog[0].size()==cd_.rhog[1].size());
  const int ng = cd_.rhog[0].size();
  vector<complex<double> > rhog_in(nspin*ng), drhog(nspin*ng),
    rhobar(nspin*ng), drhobar(nspin*ng);
  const int anderson_ndim = s_.ctrl.charge_mix_ndim;
  vector<double> wkerker(ng), wls(ng);

  // Anderson charge mixer: include both spins in the same vector
  // Factor of 2: complex coeffs stored as double
  MPI_Comm vcomm = MPIdata::g_comm();
  AndersonMixer mixer(2*nspin*ng,anderson_ndim,&vcomm);

  // compute Kerker preconditioning
  // real space Kerker cutoff in a.u.
  const double rc_Kerker = s_.ctrl.charge_mix_rcut;
  const double *const g2 = cd_.vbasis()->g2_ptr();

  // define q1 cutoff for row weighting of LS charge mixing
  // Use rc1 = 3 a.u. default cutoff
  double rc1 = 3.0;
  // check if override from the debug map
  // use: set debug RC1 <value>
  imap = debug_map.find("RC1");
  if ( imap != debug_map.end() )
  {
    const string val = imap->second;
    istringstream is(val);
    is >> rc1;
    if ( onpe0 )
      cout << " override rc1 value: rc1 = " << rc1 << endl;
    assert(rc1 >= 0.0);
  }

  double delta_ratio = 0.01;
  imap = debug_map.find("DELTA_RATIO");
  if ( imap != debug_map.end() )
  {
    const string val = imap->second;
    istringstream is(val);
    is >> delta_ratio;
    if ( onpe0 )
      cout << " override delta_ratio value = " << delta_ratio << endl;
    assert(delta_ratio >= 0.0);
  }

  if ( rc1 != 0.0 )
  {
    const double q1 = 2.0 * M_PI / rc1;
    for ( int i = 0; i < wls.size(); i++ )
    {
      if ( g2[i] != 0.0 )
        wls[i] = sqrt(g2[i] / ( g2[i] + q1*q1 ));
      else
        wls[i] = 1.0;
    }
  }
  else
  {
    for ( int i = 0; i < wls.size(); i++ )
      wls[i] = 1.0;
  }

  if ( rc_Kerker > 0.0 )
  {
    const double q0_kerker = 2 * M_PI / rc_Kerker;
    const double q0_kerker2 = q0_kerker * q0_kerker;
    for ( int i = 0; i < wkerker.size(); i++ )
    {
      if ( g2[i] != 0.0 )
        wkerker[i] = g2[i] / ( g2[i] + q0_kerker2 );
      else
        wkerker[i] = 1.0;
    }
  }
  else
  {
    for ( int i = 0; i < wkerker.size(); i++ )
      wkerker[i] = 1.0;
  }

  if ( onpe0 )
    cout << "<net_charge> " << atoms.nel()-wf.nel() << " </net_charge>\n";

  // Next line: special case of niter=0: compute GS only
  bool iter_done = false;
  for ( int iter = 0; iter < max(niter,1) && !iter_done; iter++ )
  {
    // ionic iteration

    tm_iter.start();

    if ( onpe0 )
      cout << "<iteration count=\"" << iter+1 << "\">\n";

    // compute energy and ionic forces using existing wavefunction
    double maxforce = 0.0;
    double maxstress = 0.0;

    if ( !gs_only )
    {
      tmap["charge"].start();
      cd_.update_density();
      tmap["charge"].stop();

      tmap["update_vhxc"].start();
      ef_.update_vhxc(compute_stress);
      tmap["update_vhxc"].stop();
      const bool compute_forces = true;
      tmap["energy"].start();
      double energy =
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
      tmap["energy"].stop();
      double enthalpy = ef_.enthalpy();

      if ( force_tol > 0.0 )
      {
        maxforce = 0.0;
        for ( int is = 0; is < fion.size(); is++ )
          for ( int i = 0; i < fion[is].size(); i++ )
            maxforce = max(maxforce, fabs(fion[is][i]));
      }

      if ( stress_tol > 0.0 )
      {
        compute_sigma();
        for ( int i = 0; i < sigma.size(); i++ )
          maxstress = max(maxstress, gpa*fabs(sigma[i]));
      }

      if ( onpe0 )
      {
        cout << cd_;
        cout << ef_;
        if ( ef_.el_enth() )
          cout << *ef_.el_enth();
      }

      if ( iter > 0 && ionic_stepper )
      {
        ionic_stepper->compute_v(energy,fion);
      }
      // at this point, positions r0, velocities v0 and forces fion are
      // consistent

      // execute commands in iter_cmd if defined
      if ( !iter_cmd_.empty() )
      {
        if ( iter % iter_cmd_period_ == 0 )
        {
          // copy positions and velocities from IonicStepper to AtomSet
          if ( ionic_stepper )
          {
            ionic_stepper->set_positions();
            ionic_stepper->set_velocities();
          }
          // command must be terminated with \n
          istringstream cmdstream(iter_cmd_ + "\n");
          s_.ui->processCmds(cmdstream,"[iter_cmd]",true);
          // copy positions and velocities back from AtomSet
          if ( ionic_stepper )
          {
            ionic_stepper->get_positions();
            ionic_stepper->get_velocities();
          }
        }
      }

      double ekin_ion = 0.0, temp_ion = 0.0;
      if ( ionic_stepper )
      {
        ekin_ion = ionic_stepper->ekin();
        temp_ion = ionic_stepper->temp();
      }

      // print positions, velocities and forces at time t0
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
                 << "    <position> " << pa->position() << " </position>\n"
                 << "    <velocity> " << pa->velocity() << " </velocity>\n"
                 << "    <force> "
                 << fion[is][i] << " "
                 << fion[is][i+1] << " "
                 << fion[is][i+2]
                 << " </force>\n";
            cout << "  </atom>" << endl;
            i += 3;
          }
        }
        cout << "</atomset>" << endl;
        cout << setprecision(6);
        cout << "<unit_cell_a_norm> " << atoms.cell().a_norm(0)
             << " </unit_cell_a_norm>" << endl;
        cout << "<unit_cell_b_norm> " << atoms.cell().a_norm(1)
             << " </unit_cell_b_norm>" << endl;
        cout << "<unit_cell_c_norm> " << atoms.cell().a_norm(2)
             << " </unit_cell_c_norm>" << endl;
        cout << setprecision(3) << "<unit_cell_alpha>  "
             << atoms.cell().alpha() << " </unit_cell_alpha>" << endl;
        cout << setprecision(3) << "<unit_cell_beta>   "
             << atoms.cell().beta() << " </unit_cell_beta>" << endl;
        cout << setprecision(3) << "<unit_cell_gamma>  "
             << atoms.cell().gamma() << " </unit_cell_gamma>" << endl;
        cout << setprecision(3) << "<unit_cell_volume> "
             << atoms.cell().volume() << " </unit_cell_volume>" << endl;

        // include the kinetic energy of the stepper
        // e.g. to include thermostat contributions
        double ekin_stepper;
        if ( ionic_stepper != 0 )
          ekin_stepper = ionic_stepper->ekin_stepper();
        cout << setprecision(8);
        cout << "  <econst> " << enthalpy+ekin_ion+ekin_stepper
             << " </econst>\n";
        cout << "  <ekin_ion> " << ekin_ion << " </ekin_ion>\n";
        cout << "  <temp_ion> " << temp_ion << " </temp_ion>\n";
      }

      if ( atoms_move )
      {
        if ( s_.constraints.size() > 0 )
        {
          s_.constraints.update_constraints(dt);
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
        compute_sigma();
        print_stress();

        if ( cell_moves )
        {
          cell_stepper->compute_new_cell(enthalpy,sigma,fion);

          // Update cell
          cell_stepper->update_cell();

          ef_.cell_moved();
          ef_.atoms_moved(); // modifications of the cell also move ions
        }
      }
    } // if !gs_only

    // Recalculate ground state wavefunctions
    // wavefunction extrapolation
    if ( atoms_move && extrapolate_wf )
    {
      for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
        {
          if ( ntc_extrapolation )
          {
            SlaterDet* sd = wf.sd(isp_loc,ikp_loc);
            SlaterDet* sdv = wfv->sd(isp_loc,ikp_loc);
            SlaterDet* sdmm = wfmm->sd(isp_loc,ikp_loc);
            double* c = (double*) sd->c().cvalptr();
            double* cv = (double*) sdv->c().cvalptr();
            double* cmm = (double*) sdmm->c().cvalptr();
            const int mloc = sd->c().mloc();
            const int nloc = sd->c().nloc();
            const int len = 2*mloc*nloc;
            if ( iter == 0 )
            {
              // copy c on cv
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double v = cv[i];
                // extrapolation using velocity in cv
                c[i] = x + dt * v;
                cv[i] = x;
              }
              tmap["gram"].start();
              sd->gram();
              tmap["gram"].stop();
            }
            else if ( iter == 1 )
            {
              sdv->align(*sd);
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xm = cv[i];
                c[i] = 2.0 * x - xm;
                cv[i] = x;
                cmm[i] = xm;
              }
              tmap["gram"].start();
              sd->gram();
              tmap["gram"].stop();
            }
            else
            {
              // align wf with wfmm before extrapolation
              sd->align(*sdmm);

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
              sd->gram();
              tmap["gram"].stop();

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
            SlaterDet* sd = wf.sd(isp_loc,ikp_loc);
            SlaterDet* sdv = wfv->sd(isp_loc,ikp_loc);
            SlaterDet* sdmm = wfmm->sd(isp_loc,ikp_loc);
            double* c = (double*) sd->c().cvalptr();
            double* cv = (double*) sdv->c().cvalptr();
            double* cmm = (double*) sdmm->c().cvalptr();
            const int mloc = s_.wf.sd(isp_loc,ikp_loc)->c().mloc();
            const int nloc = s_.wf.sd(isp_loc,ikp_loc)->c().nloc();
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
              sd->gram();
              tmap["gram"].stop();
            }
            else if ( iter == 1 )
            {
              sdv->align(*sd);
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xm = cv[i];
                c[i] = 2.0 * x - xm;
                cv[i] = x;
                cmm[i] = xm;
              }
              tmap["gram"].start();
              sd->gram();
              tmap["gram"].stop();
            }
            else
            {
              // align wf with wfmm before extrapolation
              sd->align(*sdmm);

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
              sd->gram();
              tmap["gram"].stop();

              // c[i] now contains the extrapolated value
            }
            // c[i] is now ready for electronic iterations
          }
          else // normal extrapolation
          {
            SlaterDet* sd = wf.sd(isp_loc,ikp_loc);
            SlaterDet* sdv = wfv->sd(isp_loc,ikp_loc);
            double* c = (double*) sd->c().cvalptr();
            double* cv = (double*) sdv->c().cvalptr();
            const int mloc = s_.wf.sd(isp_loc,ikp_loc)->c().mloc();
            const int nloc = s_.wf.sd(isp_loc,ikp_loc)->c().nloc();
            const int len = 2*mloc*nloc;
            if ( iter == 0 )
            {
              // copy c to cv
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double v = cv[i];
                c[i] = x + dt * v;
                cv[i] = x;
              }
              tmap["gram"].start();
              sd->gram();
              tmap["gram"].stop();
            }
            else
            {
              tmap["align"].start();
              sdv->align(*sd);
              tmap["align"].stop();

              // linear extrapolation
              for ( int i = 0; i < len; i++ )
              {
                const double x = c[i];
                const double xm = cv[i];
                c[i] = 2.0 * x - xm;
                cv[i] = x;
              }
              tmap["lowdin"].start();
              sd->lowdin();
              tmap["lowdin"].stop();
            }
          }
        }
      }
    } // atoms_move && extrapolate_wf

    // do nitscf self-consistent iterations, each with nite electronic steps
    if ( wf_stepper != 0 )
    {
      wf_stepper->preprocess();
      if ( anderson_charge_mixing )
        mixer.restart();

      double ehart, ehart_m;
      double delta_ehart;
      bool scf_converged = false;
      int itscf = 0;
      double etotal = 0.0, etotal_m = 0.0, etotal_mm = 0.0;

      while ( !scf_converged && itscf < nitscf_ )
      {
        if ( nite_ > 0 && onpe0 )
          cout << "  BOSampleStepper: start scf iteration" << endl;

        // update charge density
        tmap["charge"].start();
        // The density is updated at the first step if update_density_first_
        // is true.
        // It is always updated after the first step
        if ( ( update_density_first_ || itscf>0 ) )
          cd_.update_density();
        tmap["charge"].stop();

        if ( onpe0 )
          cout << cd_;

        // charge mixing
        if ( nite_ > 0 )
        {
          if ( itscf == 0 )
          {
            // at first scf iteration, nothing to mix
            // memorize rhog in rhog_in for next iteration

            for ( int ispin = 0; ispin < nspin; ispin++ )
              for ( int i = 0; i < ng; i++ )
                rhog_in[i+ng*ispin] = cd_.rhog[ispin][i];
          }
          else
          {
            // itscf > 0
            // compute unscreened correction drhog
            for ( int ispin = 0; ispin < nspin; ispin++ )
            {
              for ( int i = 0; i < ng; i++ )
              {
                drhog[i+ng*ispin] = (cd_.rhog[ispin][i]-rhog_in[i+ng*ispin]);
              }
            }

            const double alpha = s_.ctrl.charge_mix_coeff;
            // Anderson acceleration
            if ( anderson_charge_mixing )
            {
              // row weighting of LS calculation
              for ( int ispin = 0; ispin < nspin; ispin++ )
              {
                for ( int i = 0; i < ng; i++ )
                  drhog[i+ng*ispin] /= wls[i];
              }

              if ( sd_ctxt.mycol() == 0 )
              {
                // use AndersonMixer on first column only and bcast results
                mixer.update((double*)&rhog_in[0], (double*)&drhog[0],
                             (double*)&rhobar[0], (double*)&drhobar[0]);
                const int n = 2*nspin*ng;
                sd_ctxt.dbcast_send('r',n,1,(double*)&rhobar[0],n);
                sd_ctxt.dbcast_send('r',n,1,(double*)&drhobar[0],n);
              }
              else
              {
                const int n = 2*nspin*ng;
                sd_ctxt.dbcast_recv('r',n,1,(double*)&rhobar[0],n,-1,0);
                sd_ctxt.dbcast_recv('r',n,1,(double*)&drhobar[0],n,-1,0);
              }

              for ( int ispin = 0; ispin < nspin; ispin++ )
              {
                for ( int i = 0; i < ng; i++ )
                  drhobar[i+ng*ispin] *= wls[i];
                for ( int i = 0; i < ng; i++ )
                  rhog_in[i+ng*ispin] = rhobar[i+ng*ispin] +
                    alpha * drhobar[i+ng*ispin] * wkerker[i];
              }
            }
            else
            {
              for ( int ispin = 0; ispin < nspin; ispin++ )
              {
                for ( int i = 0; i < ng; i++ )
                  rhog_in[i+ng*ispin] += alpha * drhog[i+ng*ispin] * wkerker[i];
              }
            }

            for ( int ispin = 0; ispin < nspin; ispin++ )
            {
              for ( int i = 0; i < ng; i++ )
                cd_.rhog[ispin][i] = rhog_in[i+ng*ispin];
            }
            cd_.update_rhor();
          }
        } // if nite_ > 0

        // update vhxc:
        // at first scf step:
        // - update both vh and vxc
        // at later steps:
        // - update depending of values of update_vh_ and update_vxc_
        tmap["update_vhxc"].start();
        if ( itscf == 0 )
          ef_.update_vhxc(compute_stress);
        else
          ef_.update_vhxc(compute_stress, update_vh_, update_vxc_);
        tmap["update_vhxc"].stop();

        // reset stepper only if multiple non-selfconsistent steps
        if ( nite_ > 0 ) wf_stepper->preprocess();

        // non-self-consistent loop
        // repeat until the change in eigenvalue_sum is smaller than a
        // fraction delta_ratio of the change in Hartree energy delta_ehart
        // in the last scf iteration
        bool nonscf_converged = false;
        if ( itscf == 0 )
        {
          ehart = ef_.ehart();
          delta_ehart = 0.0;
        }
        else if ( itscf == 1 )
        {
          ehart_m = ehart;
          ehart = ef_.ehart();
          delta_ehart = fabs(ehart - ehart_m);
        }
        else
        {
          // itscf > 1
          // only allow decrease in delta_ehart
          ehart_m = ehart;
          ehart = ef_.ehart();
          delta_ehart = min(delta_ehart,fabs(ehart - ehart_m));
        }
#if DEBUG
        if ( onpe0 )
        {
          cout << " BOSampleStepper::step: delta_ehart: "
               << delta_ehart << endl;
        }
#endif

        int ite = 0;
        double etotal_int;

        double eigenvalue_sum, eigenvalue_sum_m = 0.0;
        // if nite == 0: do 1 iteration, no screening in charge mixing
        // if nite > 0: do nite iterations, use screening in charge mixing
        //
        while ( !nonscf_converged && ite < max(nite_,1) )
        {
          tmap["energy"].start();
          ef_.energy(true,dwf,false,fion,false,sigma_eks);
          tmap["energy"].stop();

          // compute the sum of eigenvalues (with fixed weight)
          // to measure convergence of the subspace update
          // compute trace of the Hamiltonian matrix Y^T H Y
          // scalar product of Y and (HY): tr Y^T (HY) = sum_ij Y_ij (HY)_ij
          // Note: since the hamiltonian is hermitian and dwf=H*wf
          // the dot product in the following line is real

          if ( ite > 0 )
            eigenvalue_sum_m = eigenvalue_sum;

          eigenvalue_sum = real(s_.wf.dot(dwf));
          if ( onpe0 )
            cout << "  <eigenvalue_sum>  " << setprecision(8)
                 << eigenvalue_sum << " </eigenvalue_sum>" << endl;

          tmap["wf_update"].start();
          wf_stepper->update(dwf);
          tmap["wf_update"].stop();

          // compare delta_eig_sum only after first iteration
          if ( ite > 0 )
          {
            double delta_eig_sum = fabs(eigenvalue_sum - eigenvalue_sum_m);
            nonscf_converged |= (delta_eig_sum < delta_ratio * delta_ehart);
#if DEBUG
            if ( onpe0 )
            {
              cout << " BOSampleStepper::step delta_eig_sum: "
                   << delta_eig_sum << endl;
            }
#endif
          }
          ite++;
        }

        // subspace diagonalization
        if ( compute_eigvec || s_.ctrl.wf_diag == "EIGVAL" )
        {
          tmap["energy"].start();
          ef_.energy(true,dwf,false,fion,false,sigma_eks);
          tmap["energy"].stop();
          tmap["wfdiag"].start();
          s_.wf.diag(dwf,compute_eigvec);
          tmap["wfdiag"].stop();
          // print eigenvalues
          if ( MPIdata::onpe0() )
            cout << "<eigenset>" << endl;
          for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
          {
            const int isp_loc = wf.isp_local(ispin);
            for ( int ikp = 0; ikp < wf.nkp(); ++ikp )
            {
              const int ikp_loc = wf.ikp_local(ikp);
              ostringstream ostr;
              ostr.setf(ios::fixed,ios::floatfield);
              ostr.setf(ios::right,ios::adjustfield);
              int isrc = -1;
              if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
              {
                if ( MPIdata::sd_rank() == 0 )
                {
                  ostr.str("");
                  isrc = MPIdata::rank();
                  const int nst = wf.sd(isp_loc,ikp_loc)->nst();
                  const double eVolt = 2.0 * 13.6058;
                  ostr <<    "  <eigenvalues spin=\"" << ispin
                       << "\" kpoint=\""
                       << setprecision(8)
                       << wf.sd(isp_loc,ikp_loc)->kpoint()
                       << "\" weight=\""
                       << setprecision(8)
                       << wf.weight(ikp)
                       << "\" n=\"" << nst << "\">" << endl;
                  for ( int i = 0; i < nst; i++ )
                  {
                    ostr << setw(12) << setprecision(5)
                         << wf.sd(isp_loc,ikp_loc)->eig(i)*eVolt;
                    if ( i%5 == 4 ) ostr << endl;
                  }
                  if ( nst%5 != 0 ) ostr << endl;
                  ostr << "  </eigenvalues>" << endl;
                }
              }
              cout0(ostr.str(),isrc);
              MPI_Barrier(MPIdata::comm());
            }
          }
          if ( MPIdata::onpe0() )
            cout << "</eigenset>" << endl;
          MPI_Barrier(MPIdata::comm());
        }

        // update occupation numbers if fractionally occupied states
        // compute weighted sum of eigenvalues
        // default value if no fractional occupation
        double fac = ( wf.nspin() == 1 ) ? 2.0 : 1.0 ;
        double w_eigenvalue_sum = fac * eigenvalue_sum;

        if ( fractional_occ )
        {
          if ( s_.ctrl.fermi_temp > 0 )
            wf.update_occ(s_.ctrl.fermi_temp);
#if 0
          if ( onpe0 )
          {
            cout << "  Wavefunction entropy: " << wf_entropy << endl;
            cout << "  Entropy contribution to free energy: "
                 << - wf_entropy * s_.ctrl.fermi_temp * boltz << endl;
          }
#endif
          w_eigenvalue_sum = 0.0;
          for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
          {
            for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
            {
              const int nst = wf.sd(isp_loc,ikp_loc)->nst();
              const int ikpg = wf.ikp_global(ikp_loc);
              const double wkp = wf.weight(ikpg);
              for ( int n = 0; n < nst; n++ )
              {
                const double occ = wf.sd(isp_loc,ikp_loc)->occ(n);
                w_eigenvalue_sum += wkp * occ * wf.sd(isp_loc,ikp_loc)->eig(n);
              }
            }
          }
          double tsum;
          MPI_Allreduce(&w_eigenvalue_sum,&tsum,1,
            MPI_DOUBLE,MPI_SUM,MPIdata::kp_sp_comm());
          w_eigenvalue_sum = tsum;
        }

        // Harris-Foulkes estimate of the total energy
        etotal_int = w_eigenvalue_sum - ef_.ehart_e() + ef_.ehart_p() +
                     ef_.esr() - ef_.eself() + ef_.dxc() + ef_.ets();
#ifdef DEBUG
        if ( onpe0 )
        {
          cout << setprecision(8);
          cout << "w_eigenvalue_sum = " << setw(15) << w_eigenvalue_sum << endl;
          cout << "ef.dxc()         = " << setw(15) << ef_.dxc() << endl;
          cout << "ef.ehart()       = " << setw(15) << ef_.ehart() << endl;
          cout << "ef.ehart_e()     = " << setw(15) << ef_.ehart_e() << endl;
          cout << "ef.ehart_ep()    = " << setw(15) << ef_.ehart_ep() << endl;
          cout << "ef.esr()         = " << setw(15) << ef_.esr() << endl;
        }
#endif

        if ( onpe0 )
        {
          cout.setf(ios::fixed,ios::floatfield);
          cout.setf(ios::right,ios::adjustfield);
          cout << "  <etotal_int>  " << setprecision(8) << setw(15)
               << etotal_int << " </etotal_int>\n";
        }

        etotal_mm = etotal_m;
        etotal_m = etotal;
        etotal = etotal_int;

        if ( nite_ > 0 && onpe0 )
          cout << "  BOSampleStepper: end scf iteration" << endl;

        // delta_etotal = interval containing etotal, etotal_m and etotal_mm
        double delta_etotal = fabs(etotal - etotal_m);
        delta_etotal = max(delta_etotal,fabs(etotal - etotal_mm));
        delta_etotal = max(delta_etotal,fabs(etotal_m - etotal_mm));
        // bcast the value of delta_etotal to ensure coherence
        // of scf_converged
        MPI_Bcast(&delta_etotal,1,MPI_DOUBLE,0,MPIdata::comm());
        scf_converged |= (delta_etotal < s_.ctrl.scf_tol);
        itscf++;
      } // while scf

      if ( compute_mlwf || compute_mlwfc )
      {
        tmap["mlwf"].start();
        for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
        {
          const int isp_loc = wf.isp_local(ispin);
          const int ikp = 0;
          const int ikp_loc = wf.ikp_local(ikp);
          if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
          {
            mlwft[isp_loc]->update();
            mlwft[isp_loc]->compute_transform();

            if ( compute_mlwf )
            {
              SlaterDet& sd = *(wf.sd(isp_loc,ikp_loc));
              mlwft[isp_loc]->apply_transform(sd);
            }
          }
        }

        // print mlwf centers
        D3vector edipole_sum, d3tsum;
        if ( onpe0 )
          cout << "<mlwfs>" << endl;
        for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
        {
          const int isp_loc = wf.isp_local(ispin);
          const int ikp = 0;
          const int ikp_loc = wf.ikp_local(ikp);
          ostringstream ostr;
          int isrc = -1;
          if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
          {
            SlaterDet& sd = *(wf.sd(isp_loc,ikp_loc));
            if ( MPIdata::sd_rank() == 0 )
            {
              ostr.str("");
              ostr.setf(ios::fixed,ios::floatfield);
              ostr.setf(ios::right,ios::adjustfield);
              isrc = MPIdata::rank();
              ostr << " <mlwfset spin=\"" << ispin
                   << "\" size=\"" << sd.nst() << "\">" << endl;
              double total_spread[6];
              for ( int j = 0; j < 6; j++ )
                 total_spread[j] = 0.0;
              for ( int i = 0; i < sd.nst(); i++ )
              {
                D3vector ctr = mlwft[isp_loc]->center(i);
                double spi[6];
                for (int j=0; j<3; j++)
                {
                  spi[j] = mlwft[isp_loc]->spread2(i,j);
                  total_spread[j] += spi[j];
                }

                ostr << "   <mlwf center=\"" << setprecision(6)
                     << setw(12) << ctr.x
                     << setw(12) << ctr.y
                     << setw(12) << ctr.z
                     << " \" spread=\" "
                     << setw(12) << spi[0]
                     << setw(12) << spi[1]
                     << setw(12) << spi[2] << " \"/>"
                     << endl;
              }

              ostr << " <total_spread> ";
              for ( int j = 0; j < 3; j++ )
                ostr << setprecision(6) << setw(15) << total_spread[j];
              ostr << " </total_spread>" << endl;
              D3vector edipole = mlwft[isp_loc]->dipole();
              ostr << " <electronic_dipole spin=\"" << ispin
                   << "\"> " << edipole << " </electronic_dipole>" << endl;
              edipole_sum += edipole;
              ostr << " </mlwfset>" << endl;
            } // sd_rank() == 0
          }
          cout0(ostr.str(),isrc);
          MPI_Barrier(MPIdata::comm());
        } // ispin
        if ( onpe0 )
          cout << "</mlwfs>" << endl;

        if ( MPIdata::sd_rank() == 0 )
        {
          MPI_Reduce(&edipole_sum[0],&d3tsum[0],3,
                     MPI_DOUBLE,MPI_SUM,0,MPIdata::sp_comm());
          edipole_sum = d3tsum;
        }

        if ( onpe0 )
        {
          D3vector idipole = atoms.dipole();
          cout << setprecision(6);
          cout << " <ionic_dipole> " << idipole
               << " </ionic_dipole>" << endl;
          cout << " <total_dipole> " << idipole + edipole_sum
               << " </total_dipole>" << endl;
          cout << " <total_dipole_length> " << length(idipole + edipole_sum)
               << " </total_dipole_length>" << endl;
        }
        tmap["mlwf"].stop();
      }

      // If GS calculation only, print energy and atomset at end of iterations
      if ( gs_only )
      {
        tmap["charge"].start();
        cd_.update_density();
        tmap["charge"].stop();

        tmap["update_vhxc"].start();
        ef_.update_vhxc(compute_stress);
        tmap["update_vhxc"].stop();
        const bool compute_forces = true;
        tmap["energy"].start();
        ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
        tmap["energy"].stop();

        if ( onpe0 )
        {
          cout << ef_;
          if ( ef_.el_enth() )
            cout << *ef_.el_enth();
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
                   << "    <position> " << pa->position() << " </position>\n"
                   << "    <velocity> " << pa->velocity() << " </velocity>\n"
                   << "    <force> "
                   << fion[is][i] << " "
                   << fion[is][i+1] << " "
                   << fion[is][i+2]
                   << " </force>\n";
              cout << "  </atom>" << endl;
              i += 3;
            }
          }
          cout << "</atomset>" << endl;
          cout << setprecision(6);
          cout << "<unit_cell_a_norm> " << atoms.cell().a_norm(0)
               << " </unit_cell_a_norm>" << endl;
          cout << "<unit_cell_b_norm> " << atoms.cell().a_norm(1)
               << " </unit_cell_b_norm>" << endl;
          cout << "<unit_cell_c_norm> " << atoms.cell().a_norm(2)
               << " </unit_cell_c_norm>" << endl;
          cout << setprecision(3) << "<unit_cell_alpha>  "
               << atoms.cell().alpha() << " </unit_cell_alpha>" << endl;
          cout << setprecision(3) << "<unit_cell_beta>   "
               << atoms.cell().beta() << " </unit_cell_beta>" << endl;
          cout << setprecision(3) << "<unit_cell_gamma>  "
               << atoms.cell().gamma() << " </unit_cell_gamma>" << endl;
          cout << setprecision(3) << "<unit_cell_volume> "
               << atoms.cell().volume() << " </unit_cell_volume>" << endl;
          if ( compute_stress )
          {
            compute_sigma();
            print_stress();
          }
        }
      }
      wf_stepper->postprocess();
    }
    else
    {
      // wf_stepper == 0, wf_dyn == LOCKED
      // evaluate and print energy
      tmap["charge"].start();
      cd_.update_density();
      tmap["charge"].stop();
      tmap["update_vhxc"].start();
      ef_.update_vhxc(compute_stress);
      tmap["update_vhxc"].stop();
      tmap["energy"].start();
      ef_.energy(true,dwf,false,fion,false,sigma_eks);
      tmap["energy"].stop();
      if ( onpe0 )
      {
        cout << cd_;
        cout << ef_;
        if ( ef_.el_enth() )
          cout << *ef_.el_enth();
      }
    }

    // if using force_tol or stress_tol, check if maxforce and maxstress
    // within tolerance
    if ( force_tol > 0.0 )
    {
      if ( onpe0 )
        cout << "  maxforce: " << scientific
             << setprecision(4) << maxforce << fixed << endl;
      iter_done |= ( maxforce < force_tol );
    }
    if ( stress_tol > 0.0 )
    {
      if ( onpe0 )
        cout << "  maxstress: " << scientific
             << setprecision(4) << maxstress << fixed << endl;
      iter_done |= ( maxstress < stress_tol );
    }

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

  } // for iter

  if ( atoms_move )
  {
    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    tmap["update_vhxc"].start();
    ef_.update_vhxc(compute_stress);
    tmap["update_vhxc"].stop();
    const bool compute_forces = true;
    tmap["energy"].start();
    double energy =
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
    tmap["energy"].stop();

    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
  }

  if ( atoms_move && extrapolate_wf )
  {
    // compute wavefunction velocity after last iteration
    // s_.wfv contains the previous wavefunction

    tmap["align"].start();
    s_.wfv->align(s_.wf);
    tmap["align"].stop();

    for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
      {
        double* c = (double*) s_.wf.sd(isp_loc,ikp_loc)->c().cvalptr();
        double* cm = (double*) s_.wfv->sd(isp_loc,ikp_loc)->c().cvalptr();
        const int mloc = s_.wf.sd(isp_loc,ikp_loc)->c().mloc();
        const int nloc = s_.wf.sd(isp_loc,ikp_loc)->c().nloc();
        const int len = 2*mloc*nloc;
        const double dt_inv = 1.0 / dt;
        if ( ntc_extrapolation )
        {
          double* cmm = (double*) wfmm->sd(isp_loc,ikp_loc)->c().cvalptr();
          for ( int i = 0; i < len; i++ )
          {
            const double x = c[i];
            const double xmm = cmm[i];
            cm[i] = dt_inv * ( x - xmm );
          }
          tmap["gram"].start();
          s_.wf.sd(isp_loc,ikp_loc)->gram();
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
          s_.wf.sd(isp_loc,ikp_loc)->gram();
          tmap["gram"].stop();
        }
      }
    }

    // compute ionic forces at last position to update velocities
    // consistently with last position
    tmap["charge"].start();
    cd_.update_density();
    tmap["charge"].stop();

    tmap["update_vhxc"].start();
    ef_.update_vhxc(compute_stress);
    tmap["update_vhxc"].stop();
    const bool compute_forces = true;
    tmap["energy"].start();
    double energy =
      ef_.energy(false,dwf,compute_forces,fion,compute_stress,sigma_eks);
    tmap["energy"].stop();

    ionic_stepper->compute_v(energy,fion);
    // positions r0 and velocities v0 are consistent
  }

  for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
    delete mlwft[isp_loc];

  // delete steppers
  delete wf_stepper;
  delete ionic_stepper;
  delete cell_stepper;

  if ( ntc_extrapolation || asp_extrapolation ) delete wfmm;

  update_density_first_ = true;
}
