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
// ResponseCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include<cassert>
#include<vector>
#include<string>
#include<fstream>
using namespace std;

#include "BOSampleStepper.h"
#include "ExternalPotential.h"
#include "FourierTransform.h"
#include "ResponseCmd.h"
#include "release.h"
#include "isodate.h"
#include "Species.h"

////////////////////////////////////////////////////////////////////////////////
int ResponseCmd::action(int argc, char **argv)
{
  // " syntax: response amplitude nitscf [nite]\n\n"
  // " syntax: response -vext vext_file [-RPA] [-amplitude a] nitscf [nite]\n\n"
  if ( argc < 3 || argc > 8)
  {
    if ( ui->onpe0() )
      cout << "  use: response amplitude nitscf [nite]\n"
           << "       response -vext vext_file [-RPA]"
           << " [-amplitude a] nitscf [nite]\n"
           << endl;
    return 1;
  }
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

  int nitscf;
  int nite = 0;
  bool rpa = false;
  int iarg = 1;
  if ( !strcmp(argv[iarg],"-vext") )
  {
    iarg++;
    // response to vext
    if ( s->vext )
      delete s->vext;
    s->vext = new ExternalPotential(*s,argv[iarg]);
    iarg++;

    if ( !strcmp(argv[iarg],"-RPA") )
    {
      rpa = true;
      iarg++;
    }

    if ( !strcmp(argv[iarg],"-amplitude") )
    {
      iarg++;
      double amplitude = atof(argv[iarg]);
      s->vext->set_amplitude(amplitude);
      iarg++;
    }

    nitscf = atoi(argv[iarg]);
    iarg++;

    if ( iarg < argc )
      nite = atoi(argv[iarg]);

    if ( nitscf == 0 )
    {
      if ( ui->onpe0() )
        cout << "  ResponseCmd: nitscf = 0, cannot run" << endl;
      return 1;
    }

    responseVext(rpa,nitscf,nite);
  }
  else
  {
    // polarizability calculation
    assert(argc > 2);
    if ( ui->onpe0() )
      cout << " ResponseCmd: polarizability with "
           << " amplitude " << argv[1] << endl;
    double amplitude = atof(argv[1]);
    nitscf = atoi(argv[2]);
    if ( argc == 4 )
      nite = atoi(argv[3]);
    responseEfield(amplitude,nitscf,nite);
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void ResponseCmd::responseEfield(double amplitude, int nitscf, int nite)
{
  // compute dipole change
  SampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);
  assert(stepper!=0);
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
void ResponseCmd::responseVext(bool rpa, int nitscf, int nite)
{
  const int nspin = s->wf.nspin();
  BOSampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);

  if ( rpa )
    stepper->set_update_vxc(false);

  // save a copy of initial wave functions
  Wavefunction wf0(s->wf);
  wf0 = s->wf;

  ChargeDensity &cd = stepper->cd();
  MPI_Comm vcomm = cd.vcomm();
  int mype, vcomm_size;
  MPI_Comm_rank(vcomm,&mype);
  MPI_Comm_size(vcomm,&vcomm_size);

  stepper->step(0);
  const vector<vector<double> > rhor1 = cd.rhor;  // density with +Vext

  s->wf = wf0;
  s->vext->reverse();

  stepper->step(0);
  const vector<vector<double> > rhor2 = cd.rhor;  // density with -Vext
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
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    drho_r[ispin].resize(np012loc);
    drho_r_tmp[ispin].resize(np012loc);
    for ( int ir = 0; ir < np012loc; ir++ )
    {
      drho_r[ispin][ir] = rhor1[ispin][ir] - rhor2[ispin][ir];
      drho_r_tmp[ispin][ir] = complex<double>( drho_r[ispin][ir] * omega, 0.0);
    }
  }

  // Fourier (forward) transform drho_r to the basis compatible
  // with the vext grid size
  const UnitCell& cell = cd.vbasis()->cell();
  const double ecut = s->vext->ecut();
  const int np0 = cd.vft()->np0();
  const int np1 = cd.vft()->np1();
  const int np2 = cd.vft()->np2();
  Basis basis(vcomm,D3vector(0,0,0));
  basis.resize(cell,cell,ecut);
  FourierTransform ft1 (basis,np0,np1,np2);

  vector<vector<complex<double> > > drho_g;
  drho_g.resize(nspin);
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    drho_g[ispin].resize(basis.localsize());
    ft1.forward(&drho_r_tmp[ispin][0],&drho_g[ispin][0]);
  }

  // Fourier (backward) transform drho_g to drho_r on the
  // vext grid (n0,n1,n2)
  const int n0 = s->vext->n(0);
  const int n1 = s->vext->n(1);
  const int n2 = s->vext->n(2);
  const int n012 = n0 * n1 * n2;
  FourierTransform ft2(basis, n0, n1, n2);
  const int n012loc = ft2.np012loc();

  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    drho_r[ispin].resize(n012loc);
    drho_r_tmp[ispin].resize(n012loc);
    ft2.backward(&drho_g[ispin][0],&drho_r_tmp[ispin][0]);
    for ( int ir = 0; ir < drho_r[ispin].size(); ir++ )
    {
      drho_r[ispin][ir] = 0.5 * real( drho_r_tmp[ispin][ir] * omega_inv ) /
                          s->vext->amplitude();
    }
  }

  // collect parts of drho_r to first task on column context and write
  vector<double> rbuf;
  if ( mype == 0 )
    rbuf.resize(ft2.np012());

  vector<int> rcounts(vcomm_size,0);
  vector<int> displs(vcomm_size,0);
  int displ = 0;
  for ( int iproc = 0; iproc < vcomm_size; iproc++ )
  {
    displs[iproc] = displ;
    displ += ft2.np012loc(iproc);
    rcounts[iproc] = ft2.np012loc(iproc);
  }

  ofstream os;
  for ( int ispin=0; ispin < nspin; ispin++ )
  {
    if ( ui->onpe0() )
    {
      string filename;
      if ( nspin == 1 )
        filename = s->vext->filename() + ".response";
      else
        filename = s->vext->filename() + ".response."
                   + ( (ispin==0) ? "spin0" : "spin1" );
      os.open(filename.c_str());
    }

    MPI_Gatherv(&drho_r[ispin][0],ft2.np012loc(),MPI_DOUBLE,&rbuf[0],
                &rcounts[0],&displs[0],MPI_DOUBLE,0,vcomm);

    if  ( ui->onpe0() )
    {
      // write cube file
      os << "Created " << isodate() << " by qbox-" << release() << endl;
      os << "Charge density response under external potential "
         << s->vext->filename() << endl;
      // atoms and unit cell
      int natoms = s->atoms.size();
      D3vector a0 = s->atoms.cell().a(0);
      D3vector a1 = s->atoms.cell().a(1);
      D3vector a2 = s->atoms.cell().a(2);
      os << natoms << " " << -0.5*(a0+a1+a2) << endl;
      os << n0 << " " << a0/n0 << endl;
      os << n1 << " " << a1/n1 << endl;
      os << n2 << " " << a2/n2 << endl;
      const int nsp = s->atoms.nsp();
      for ( int is = 0; is < nsp; is++ )
      {
        Species* sp = s->atoms.species_list[is];
        const int z = sp->atomic_number();
        const int na = s->atoms.na(is);
        for ( int ia = 0; ia < na; ia++ )
        {
          Atom *ap = s->atoms.atom_list[is][ia];
          os << setprecision(5);
          os << z << " " << ((double) z) << " " << ap->position() << endl;
        }
      }
      // charge density response
      for ( int nx = 0; nx < n0; nx++ )
        for ( int ny = 0; ny < n1; ny++ )
          for ( int nz = 0; nz < n2; nz++ )
          {
            const int ir = nx + ny * n0 + nz * n0 * n1;
            os << setprecision(12) << std::scientific << rbuf[ir] << endl;
          }
      os.close();
    }
  }

  delete stepper;
}
