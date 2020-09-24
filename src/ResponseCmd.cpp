////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2016 The Regents of the University of California
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
// ResponseCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "BOSampleStepper.h"
#include "ExternalPotential.h"
#include "FourierTransform.h"
#include "ResponseCmd.h"
#include "release.h"
#include "isodate.h"
#include "Species.h"
#include "Function3d.h"
#include<cassert>
#include<vector>
#include<string>
#include<fstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int ResponseCmd::action(int argc, char **argv)
{
  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << "  ResponseCmd: no states, cannot run" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << "  ResponseCmd: ecut = 0.0, cannot run" << endl;
    return 1;
  }
  if ( s->wf.nel() != 2 * s->wf.nst() )
  {
    if ( ui->onpe0() )
      cout << "  ResponseCmd: cannot run with fractionally\n"
              "  occupied or empty states" << endl;
    return 1;
  }

  bool rpa = false;
  bool ipa = false;
  double amplitude = 0.0;

  int iarg = 1;
  if ( !strcmp(argv[iarg],"-vext") )
  {
    // response to vext
    if ( s->vext )
    {
      if ( ui->onpe0() )
        cout << "  ResponseCmd: cannot run when vext is already set" << endl;
      return 1;
    }
    iarg++;
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    string filename = argv[iarg];
    iarg++;
    string fmt = "xml";
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    if ( !strcmp(argv[iarg],"-cube") )
    {
      fmt = "cube";
      iarg++;
    }
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    if ( !strcmp(argv[iarg],"-RPA") )
    {
      rpa = true;
      iarg++;
    }
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    if ( !strcmp(argv[iarg],"-IPA") )
    {
      ipa = true;
      iarg++;
    }
    if ( rpa && ipa )
    {
      if ( ui->onpe0() )
        cout << " Only one of -RPA or -IPA can be specified" << endl;
      return 1;
    }
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    if ( !strcmp(argv[iarg],"-amplitude") )
    {
      iarg++;
      if ( iarg >= argc )
      {
        if ( ui->onpe0() )
          cout << help_msg();
        return 1;
      }
      amplitude = atof(argv[iarg]);
      iarg++;
    }

    int nitscf = 0;
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    nitscf = atoi(argv[iarg]);
    iarg++;
    int nite = 0;
    if ( iarg < argc )
      nite = atoi(argv[iarg]);

    if ( nitscf == 0 )
    {
      if ( ui->onpe0() )
        cout << "  ResponseCmd: nitscf = 0, cannot run" << endl;
      return 1;
    }

    s->vext = new ExternalPotential(*s, filename, fmt);
    s->vext->set_amplitude(amplitude);
    responseVext(rpa, ipa, nitscf, nite, fmt);
    delete s->vext;
    s->vext = 0;
  }
  else
  {
    // polarizability calculation in constant field
    // response amplitude [-RPA|-IPA] nitscf [nite]
    if ( argc < 3 || argc > 5 )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    amplitude = atof(argv[iarg]);
    iarg++;
    if ( !strcmp(argv[iarg],"-RPA") )
    {
      rpa = true;
      iarg++;
    }
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    if ( !strcmp(argv[iarg],"-IPA") )
    {
      ipa = true;
      iarg++;
    }
    if ( rpa && ipa )
    {
      if ( ui->onpe0() )
        cout << " Only one of -RPA or -IPA can be specified" << endl;
      return 1;
    }
    if ( iarg >= argc )
    {
      if ( ui->onpe0() )
        cout << help_msg();
      return 1;
    }
    int nitscf = atoi(argv[iarg]);
    iarg++;
    int nite = 0;
    if ( iarg < argc )
      nite = atoi(argv[iarg]);

    if ( ui->onpe0() )
      cout << " ResponseCmd: polarizability with "
           << " amplitude " << argv[1] << endl;
    responseEfield(amplitude,rpa,ipa,nitscf,nite);
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void ResponseCmd::responseEfield(double amplitude, bool rpa, bool ipa,
  int nitscf, int nite)
{
  // compute dipole change
  BOSampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);
  assert(stepper!=0);

  if ( rpa )
    stepper->set_update_vxc(false);

  if ( ipa )
  {
    stepper->set_update_vh(false);
    stepper->set_update_vxc(false);
  }

  ElectricEnthalpy* el_enth = stepper->ef().el_enth();
  if ( el_enth == 0 )
  {
    if ( ui->onpe0() )
      cout << "ResponseCmd:responseEfield: ElectricEnthalpy not defined\n"
           << "The polarization variable may be set to OFF" << endl;
    return;
  }

  D3vector e_field[3] = { D3vector(amplitude,0.0,0.0),
                          D3vector(0.0,amplitude,0.0),
                          D3vector(0.0,0.0,amplitude) };
  D3vector dipole_p[3], dipole_m[3];

  // save wave functions
  Wavefunction wf0(s->wf);
  wf0 = s->wf;

  // compute change in dipole in 3 directions by finite difference
  D3vector e_field_base = el_enth->e_field();
  for ( int idir = 0; idir < 3; idir++ )
  {
    el_enth->set_e_field(e_field_base+e_field[idir]);
    stepper->step(0);
    dipole_p[idir] = el_enth->dipole_total();
    s->wf = wf0;

    el_enth->set_e_field(e_field_base-e_field[idir]);
    stepper->step(0);
    dipole_m[idir] = el_enth->dipole_total();
    s->wf = wf0;
  }

  D3vector ddx = dipole_p[0] - dipole_m[0];
  D3vector ddy = dipole_p[1] - dipole_m[1];
  D3vector ddz = dipole_p[2] - dipole_m[2];

  const UnitCell& cell = s->wf.cell();
  cell.fold_in_ws(ddx);
  cell.fold_in_ws(ddy);
  cell.fold_in_ws(ddz);

  ddx *= 0.5 / amplitude;
  ddy *= 0.5 / amplitude;
  ddz *= 0.5 / amplitude;

  if ( ui->onpe0() )
  {
    cout << "<polarizability>" << endl;
    cout << setprecision(6);
    cout << " <a_xx> " << setw(10) << ddx.x << " </a_xx>"
            " <a_yx> " << setw(10) << ddx.y << " </a_yx>"
            " <a_zx> " << setw(10) << ddx.z << " </a_zx>" << endl;
    cout << " <a_xy> " << setw(10) << ddy.x << " </a_xy>"
            " <a_yy> " << setw(10) << ddy.y << " </a_yy>"
            " <a_zy> " << setw(10) << ddy.z << " </a_zy>" << endl;
    cout << " <a_xz> " << setw(10) << ddz.x << " </a_xz>"
            " <a_yz> " << setw(10) << ddz.y << " </a_yz>"
            " <a_zz> " << setw(10) << ddz.z << " </a_zz>" << endl;
    double a_iso = (ddx.x+ddy.y+ddz.z)/3.0;
    cout << " <a_iso> " << setw(10) << a_iso << " </a_iso>" << endl;
    double beta[9] = { ddx.x-a_iso, ddx.y, ddx.z,
                       ddy.x, ddy.y-a_iso, ddy.z,
                       ddz.x, ddz.y, ddz.z-a_iso };
    double a_aniso = 0.0;
    for ( int i = 0; i < 9; i++ )
      a_aniso += beta[i] * beta[i];
    a_aniso *= 2.0/15.0;
    cout << " <a_aniso> " << setw(10) << a_aniso << " </a_aniso>" << endl;
    cout << "</polarizability>" << endl;
  }

  delete stepper;
}

////////////////////////////////////////////////////////////////////////////////
void ResponseCmd::responseVext(bool rpa, bool ipa, int nitscf, int nite,
   string fmt)
{
  s->wf.info(cout, "wavefunction");

  const int nspin = s->wf.nspin();
  BOSampleStepper *stepper = new BOSampleStepper(*s, nitscf, nite);

  if ( rpa )
    stepper->set_update_vxc(false);

  if ( ipa )
  {
    stepper->set_update_vh(false);
    stepper->set_update_vxc(false);
  }

  assert(fmt == "xml" || fmt == "cube");

  // save a copy of initial wave functions
  Wavefunction wf0(s->wf);
  wf0 = s->wf;

  ChargeDensity &cd = stepper->cd();
  MPI_Comm vcomm = MPIdata::g_comm();

  stepper->step(0);
  vector<vector<double> > rhor1; // density with +Vext
  rhor1 = cd.rhor;

  s->wf = wf0;
  s->vext->reverse();

  stepper->step(0);
  vector<vector<double> > rhor2; // density with -Vext
  rhor2 = cd.rhor;
  s->vext->reverse();

  // restore initial wave functions
  s->wf = wf0;

  // compute drho_r = rhor1 - rhor2
  const int np012loc = cd.vft()->np012loc();
  const double omega = cd.vbasis()->cell().volume();
  const double omega_inv = 1.0 / omega;
  vector<vector<double> > drho_r;
  vector<vector<complex<double> > > drho_r_tmp;
  drho_r.resize(nspin);
  drho_r_tmp.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++)
  {
    drho_r[ispin].resize(np012loc);
    drho_r_tmp[ispin].resize(np012loc);
    for (int ir = 0; ir < np012loc; ir++)
    {
      drho_r[ispin][ir] = rhor1[ispin][ir] - rhor2[ispin][ir];
      drho_r_tmp[ispin][ir] = complex<double>(drho_r[ispin][ir] * omega, 0.0);
    }
  }

  // Fourier (forward) transform drho_r to the basis compatible
  // with the vext grid size
  const UnitCell &cell = cd.vbasis()->cell();
  const double ecut = s->vext->ecut();
  const int np0 = cd.vft()->np0();
  const int np1 = cd.vft()->np1();
  const int np2 = cd.vft()->np2();
  Basis basis(vcomm, D3vector(0, 0, 0));
  basis.resize(cell, cell, ecut);
  FourierTransform ft1(basis, np0, np1, np2);

  vector<vector<complex<double> > > drho_g;
  drho_g.resize(nspin);
  for (int ispin = 0; ispin < nspin; ispin++)
  {
    drho_g[ispin].resize(basis.localsize());
    ft1.forward(&drho_r_tmp[ispin][0], &drho_g[ispin][0]);
  }

  // Fourier (backward) transform drho_g to drho_r on the
  // vext grid (n0,n1,n2)
  const int n0 = s->vext->n(0);
  const int n1 = s->vext->n(1);
  const int n2 = s->vext->n(2);
  FourierTransform ft2(basis, n0, n1, n2);
  const int n012loc = ft2.np012loc();

  for (int ispin = 0; ispin < nspin; ispin++)
  {
    drho_r[ispin].resize(n012loc);
    drho_r_tmp[ispin].resize(n012loc);
    ft2.backward(&drho_g[ispin][0], &drho_r_tmp[ispin][0]);
    for (int ir = 0; ir < drho_r[ispin].size(); ir++)
    {
      drho_r[ispin][ir] = 0.5 * real(drho_r_tmp[ispin][ir] * omega_inv) /
                          s->vext->amplitude();
    }
  }
  // at this point, drho_r has been computed and distributed on column context

  // now write drho_r to disk
  const Context& ctxt = s->wf.sd_context();
  int nprow = ctxt.nprow();
  int myrow = ctxt.myrow();
  int mycol = ctxt.mycol();
  Timer tm_write_drho;

  // serial write xml file or cube file
  // processors on row 1 collect drho from other rows
  // processor at row 1, column 1 write drho to disk
  vector<double> drho_r_gathered;
  if (myrow == 0)
    drho_r_gathered.resize(ft2.np012());

  vector<int> rcounts(nprow, 0);
  vector<int> displs(nprow, 0);
  int displ = 0;
  for (int iproc = 0; iproc < nprow; iproc++)
  {
    displs[iproc] = displ;
    displ += ft2.np012loc(iproc);
    rcounts[iproc] = ft2.np012loc(iproc);
  }

  for (int ispin = 0; ispin < nspin; ispin++)
  {
    MPI_Gatherv(&drho_r[ispin][0], ft2.np012loc(), MPI_DOUBLE,
                &drho_r_gathered[0], &rcounts[0], &displs[0],
                MPI_DOUBLE, 0, vcomm);

    if ( myrow == 0 && mycol == 0 )
    {
      string filename;
      if (nspin == 1)
        filename = s->vext->filename() + ".response";
      else
        filename = s->vext->filename() + ".response."
                   + ((ispin == 0) ? "spin0" : "spin1");
      if (fmt == "xml")
      {
        tm_write_drho.start();
        Function3d f;
        f.name = "delta_rho";
        f.nx = n0;
        f.ny = n1;
        f.nz = n2;
        f.a = s->atoms.cell().a(0);
        f.b = s->atoms.cell().a(1);
        f.c = s->atoms.cell().a(2);
        f.val = drho_r_gathered;
        f.write(filename);
        tm_write_drho.stop();
      }
      else
      {
        // write cube file
        tm_write_drho.start();

        ofstream os(filename.c_str(), ios::out);
        // comment lines (first 2 lines of cube file)
        os << "Created " << isodate() << " by qbox-" << release() << "\n";
        os << "Charge density response under external potential "
           << s->vext->filename() << "\n";
        // atoms and unit cell
        int natoms = s->atoms.size();
        D3vector a0 = s->atoms.cell().a(0);
        D3vector a1 = s->atoms.cell().a(1);
        D3vector a2 = s->atoms.cell().a(2);
        os << natoms << " " << -0.5 * (a0 + a1 + a2) << "\n";
        os << n0 << " " << a0 / n0 << endl;
        os << n1 << " " << a1 / n1 << endl;
        os << n2 << " " << a2 / n2 << endl;
        const int nsp = s->atoms.nsp();
        for (int is = 0; is < nsp; is++)
        {
          Species *sp = s->atoms.species_list[is];
          const int z = sp->atomic_number();
          const int na = s->atoms.na(is);
          for (int ia = 0; ia < na; ia++)
          {
            Atom *ap = s->atoms.atom_list[is][ia];
            os << setprecision(5);
            os << z << " " << ((double) z) << " " << ap->position() << "\n";
          }
        }
        // charge density response
        os << setprecision(6) << std::scientific;
        for (int nx = 0; nx < n0; nx++)
          for (int ny = 0; ny < n1; ny++)
            for (int nz = 0; nz < n2; nz++)
            {
              const int ir = nx + ny * n0 + nz * n0 * n1;
              os << drho_r_gathered[ir] << "\n";
            }
        os.close();

        tm_write_drho.stop();
      } // if fmt
    } //if ( myrow == 0 && mycol == 0 )
  } // for ispin

  delete stepper;
}
