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
#include "mpi.h"

void writeRealSpaceFunction(const MPI_Comm&, vector<double>&, const FourierTransform&,
                            const string, const string, const bool);

////////////////////////////////////////////////////////////////////////////////
int ResponseCmd::action(int argc, char **argv)
{
  // " syntax: response amplitude nitscf [nite]\n\n"
  // " syntax: response -vext vext_file nitscf [nite]\n\n"
  if ( argc < 3 || argc > 9)
  {
    if ( ui->onpe0() )
      cout << "  use: response amplitude nitscf [nite]\n"
           << "       response -vext vext_file [-RPA or -IPA] [-amplitude amplitude] [-parallel_write] nitscf [nite]\n"
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
  int iarg = 1;
  if ( !strcmp(argv[iarg],"-vext") )
  {
    // response to vext
    Cmd *setcmd = s->ui->findCmd("set");
    iarg++; //=2
    char** setcmd_argv = new char*[2];
    setcmd_argv[0] = "set";
    setcmd_argv[1] = "vext";
    setcmd_argv[2] =  argv[iarg];
    MPI_Barrier(MPI_COMM_WORLD);
    setcmd->action(3, setcmd_argv);
    MPI_Barrier(MPI_COMM_WORLD);
    if ( ui->onpe0() )
      cout << "  ResponseCmd: compute charge density response under"
           << "  the external potential defined in " << argv[2] << endl;

    iarg++; //=3
    if ( !strcmp(argv[iarg],"-RPA") )
    {
      s->ctrl.freeze_vh = false;
      s->ctrl.freeze_vxc = true;
      if ( ui->onpe0() )
        cout << "  RPA: XC potential will be freezed"
             << "  to the value on initial density" << endl;
    }
    else if ( !strcmp(argv[iarg],"-IPA") )
    {
      s->ctrl.freeze_vh = true;
      s->ctrl.freeze_vxc = true;
      if ( ui->onpe0() )
        cout << "  IPA: Hartree and XC potential will be freezed"
             << "  to the value on initial density" << endl;
    }
    else
    {
      s->ctrl.freeze_vh = false;
      s->ctrl.freeze_vxc = false;
      if ( ui->onpe0() )
        cout << "  Hartree and XC potential will be updated" << endl;
      iarg--;
    }

    iarg++; //=3 or 4
    if ( !strcmp(argv[iarg],"-amplitude") )
    {
      iarg++; //=4 or 5
      double amplitude = atof(argv[iarg]);
      s->vext->set_amplitude(amplitude);
      if ( ui->onpe0() )
        cout << "  External potential amplitude: "
             << s->vext->amplitude() << endl;
    }
    else iarg--;

    iarg++;
    bool parallel_write = false;
    if ( !strcmp(argv[iarg],"-parallel_write") )
    {
      parallel_write = true;
      if ( ui->onpe0() )
        cout << "  Density response will be written with base64 encoding " << endl;
    }
    else iarg--;

    iarg++;
    nitscf = atoi(argv[iarg]); //iarg = 3 or 4 or 5 or 6 or 7
    //if ( ui->onpe0() )
    //  cout << " nitscf " << nitscf << endl;

    iarg++;
    if ( iarg < argc )
      nite = atoi(argv[iarg]); //iarg = 4 or 5 or 6 or 7 or 7
    //if ( ui->onpe0() )
    //  cout << " nite " << nite << endl;

    responseVext(nitscf,nite,parallel_write);

    // after response calculation, reset vext and activate hxc
    if ( ui->onpe0() )
      cout << "  End response calculation. vext is set to zero." << endl;
    setcmd_argv[2] = "NULL";
    MPI_Barrier(MPI_COMM_WORLD);
    setcmd->action(3, setcmd_argv);
    MPI_Barrier(MPI_COMM_WORLD);
    if ( s->ctrl.freeze_vh )
    {
      s->ctrl.freeze_vh = false;
      if ( ui->onpe0() )
        cout << "  Vhartree is activated again." << endl;
    }
    if ( s->ctrl.freeze_vxc )
    {
      s->ctrl.freeze_vxc = false;
      if ( ui->onpe0() )
        cout << "  Vxc is activated again." << endl;
    }
    //delete[] setcmd;
  }
  else
  {
    // polarizability calculation
    if ( ui->onpe0() )
      cout << "  ResponseCmd: start polarizability calculation with "
           << " amplitude " << argv[1] << endl;
    double amplitude = atof(argv[1]);
    nitscf = atoi(argv[2]);
    if ( argc == 3 )
      nite = atoi(argv[3]);
    cout << amplitude << nitscf << nite;
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

  D3vector e_field[3] = { D3vector(amplitude,0.0,0.0),
                          D3vector(0.0,amplitude,0.0),
                          D3vector(0.0,0.0,amplitude) };
  D3vector dipole_p[3], dipole_m[3];

  // compute change in dipole in 3 directions by finite difference
  D3vector e_field_base = el_enth->e_field();
  for ( int idir = 0; idir < 3; idir++ )
  {
    el_enth->set_e_field(e_field_base+e_field[idir]);
    stepper->step(0);
    dipole_p[idir] = el_enth->dipole_total();

    el_enth->set_e_field(e_field_base-e_field[idir]);
    stepper->step(0);
    dipole_m[idir] = el_enth->dipole_total();
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
void ResponseCmd::responseVext(int nitscf, int nite, bool parallel_write)
{
  BOSampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);
  ChargeDensity &cd = stepper->cd();
  cd.update_density();
  const vector<vector<double> > rhor0 = cd.rhor;  // GS density

  stepper->step(0);
  const vector<vector<double> > rhor1 = cd.rhor;  // density under Vext

  s->vext->reverse();
  if ( ui->onpe0() )
    cout << "  ResponseCmd: external potential is reversed,"
         << " starting another scf iteration" << endl
         << "  using 2*rho0 - rho(vext) as initial guess for charge density" << endl;

  const int nspin = s->wf.nspin();
  vector<vector<double> > rhor2_guess;
  rhor2_guess.resize(nspin);
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    rhor2_guess[ispin].resize( rhor0[ispin].size() );
    for ( int ir = 0; ir < rhor0[ispin].size(); ir++ )
    {
      rhor2_guess[ispin][ir] = 2 * rhor0[ispin][ir] - rhor1[ispin][ir];
    }
  }
  stepper->initialize_density(rhor2_guess);

  stepper->step(0);
  s->vext->reverse();
  const vector<vector<double> > rhor2 = cd.rhor;  // density under -Vext

  // compute drho_r as rhor1 - rhor2
  const int np012loc = cd.vft()->np012loc();
  const double omega = cd.vbasis()->cell().volume();
  const double omega_inv = 1.0 / omega;
  vector<vector<double> > drho_r;
  vector<vector<complex<double> > > drho_r_tmp;  // to be used for FT
  drho_r.resize(nspin);
  drho_r_tmp.resize(nspin);
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    assert(rhor1[ispin].size() == np012loc);
    assert(rhor2[ispin].size() == np012loc);
    drho_r[ispin].resize(np012loc);
    drho_r_tmp[ispin].resize(np012loc);
    for ( int ir = 0; ir < np012loc; ir++ )
    {
      drho_r[ispin][ir] = rhor1[ispin][ir] - rhor2[ispin][ir];
      drho_r_tmp[ispin][ir] = complex<double>( drho_r[ispin][ir] * omega, 0.0);
    }
  }

  // Fourier (forward) transform drho_r to the basis that is compatible
  // with the grid size of the external potential.
  const UnitCell& cell = cd.vbasis()->cell();
  const double ecut = s->vext->ecut();
  const int np0 = cd.vft()->np0();
  const int np1 = cd.vft()->np1();
  const int np2 = cd.vft()->np2();
  Basis basis(cd.vcomm(),D3vector(0,0,0));
  basis.resize(cell,cell,ecut);
  FourierTransform ft1 (basis,np0,np1,np2);

  vector<vector<complex<double> > > drho_g;
  drho_g.resize(nspin);
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    drho_g[ispin].resize(basis.localsize());
    ft1.forward(&drho_r_tmp[ispin][0],&drho_g[ispin][0]);
  }

  // Fourier (backward) transform drho_g to drho_r the grid
  // (n0,n1,n2), which is the grid of vext
  const int n0 = s->vext->n(0);
  const int n1 = s->vext->n(1);
  const int n2 = s->vext->n(2);
  const int n012 = n0 * n1 * n2;
  FourierTransform ft2 (basis, n0, n1, n2);
  const int n012loc = ft2.np012loc();

  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    drho_r[ispin].resize(n012loc);
    drho_r_tmp[ispin].resize(n012loc);
    ft2.backward(&drho_g[ispin][0],&drho_r_tmp[ispin][0]);
    for ( int ir = 0; ir < drho_r[ispin].size(); ir++ )
    {
      // return drhor scaled by 1/(2*amplitude)
      drho_r[ispin][ir] = 0.5 * real( drho_r_tmp[ispin][ir] * omega_inv ) / abs(s->vext->amplitude());
    }
  }

  const Context& ctxt = *(s->wf.spincontext());
  ctxt.barrier();
  if ( ui->onpe0() )
    cout << "  ResponseCmd: density response has been"
         << " interpolated from grid size \n"
         << " (" << np0 << ", " << np1 << ", " << np2 << ")"
         << " to (" << n0 << ", " << n1 << ", " << n2 << ") \n";

  // trace out spin degree of freedom. later we may consider spin polarized situation
  if ( nspin == 2 )
  {
    for ( int ir = 0; ir < drho_r[0].size(); ir++)
      drho_r[0][ir] += drho_r[1][ir];
    drho_r.resize(1);
  }

  if ( ctxt.mycol() == 0 )
  {
    MPI_Comm col_comm = basis.comm();
    string label = "density_response";
    string filename = s->vext->filename() + ".response";
    writeRealSpaceFunction(col_comm,drho_r[0],ft2,label,filename,parallel_write);
  }

  delete stepper;
}

void writeRealSpaceFunction(const MPI_Comm& comm, vector<double>& fr,
                            const FourierTransform& ft,
                            const string label,
                            const string filename, const bool parallel)
{
  // This function write a real space function fr that is
  // distributed among first column of processes to an xml file
  // comm is the communicator for the FIRST COLUMN of the process grid
  // ft is used to determine how fr is stored at each process
  // If parallel is true, file will be written in parallel with base64 encoding,
  // otherwise file will be written by process 0 in text format.

  assert ( fr.size() == ft.np012loc() );
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  bool onpe0 = rank == 0;

  if ( parallel )
  {
    // todo: implement MPI write, similar to SlaterDet::write()
  }
  else
  {
    // process 0 collects fr from all processes on column 0 and then write to file
    void* sbuf;
    void* rbuf;
    vector<int> rcounts(size,0);
    vector<int> displs(size,0);
    int displ = 0;
    for ( int iproc = 0; iproc < size; iproc++ )
    {
      displs[iproc] = displ;
      displ += ft.np012loc(iproc);
      rcounts[iproc] = ft.np012loc(iproc);
    }
    if ( onpe0 )
    {
      fr.resize(ft.np012(), 0.0);
      sbuf = MPI_IN_PLACE;
      rbuf = &fr[0];
    }
    else
      sbuf = &fr[0];

    MPI_Gatherv(sbuf,ft.np012loc(),MPI_DOUBLE,rbuf,
                &rcounts[0],&displs[0],MPI_DOUBLE,0,comm);
  } // if parallel


  // write output file
  stringstream header, trailer;
  header << "<" << label << " type=\"double\" nx=\"" << ft.np0()
         << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\"";
  trailer << "</" << label << ">";


  if ( parallel )
  {
    header << " encoding=\"base64\">\n";

  }
  else
  {
    if  (onpe0 )
    {
      ofstream os;
      os.open(filename.c_str());
      header << " encoding=\"text\">\n";
      os << header.str();

      for ( int ir = 0; ir < fr.size(); ir++ )
        os << setprecision(12) << std::scientific << fr[ir] << endl;

      os << trailer.str();
      os.close();
    }
  } // if parallel
  if ( onpe0 )
  {
    cout << " ResponseCmd: charge density response has been written in: "
         << filename << endl;
  }
}
