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
#include "Base64Transcoder.h"
#include <unistd.h>


////////////////////////////////////////////////////////////////////////////////
int ResponseCmd::action(int argc, char **argv)
{
  // " syntax: response amplitude nitscf [nite]\n\n"
  // " syntax: response -vext vext_file [-RPA|-IPA] [-amplitude a]
  //                    [-io iomode -nx nx -ny ny -nz nz] [-q qx qy qz] nitscf [nite]\n\n"

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
    if ( s->vext )
    {
      if ( ui->onpe0() )
        cout << "  ResponseCmd: cannot run when vext is already set" << endl;
      return 1;
    }

    string filename;
    bool rpa = false;
    bool ipa = false;
    double amplitude = 0.0;
    string io = "cube";
    int nx, ny, nz;
    nx = ny = nz = 0;

    iarg++;
    filename = argv[iarg];
    iarg++;

    if ( !strcmp(argv[iarg],"-RPA") )
    {
      rpa = true;
      iarg++;
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
    if ( !strcmp(argv[iarg],"-amplitude") )
    {
      iarg++;
      amplitude = atof(argv[iarg]);
      iarg++;
    }

    if ( !strcmp(argv[iarg],"-io") )
    {
      iarg++;
      io = argv[iarg];
      assert( io == "cube" || io == "base64_serial" || io == "base64_parallel");
      iarg++;
      if ( io == "base64_serial" || io == "base64_parallel" )
      {
        assert(!strcmp(argv[iarg],"-nx"));
        iarg++;
        nx = atoi(argv[iarg]);
        iarg++;
        assert(!strcmp(argv[iarg],"-ny"));
        iarg++;
        ny = atoi(argv[iarg]);
        iarg++;
        assert(!strcmp(argv[iarg],"-nz"));
        iarg++;
        nz = atoi(argv[iarg]);
        iarg++;
      }
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

    s->vext = new ExternalPotential(*s, filename, io, nx, ny, nz);
    s->vext->set_amplitude(amplitude);
    responseVext(rpa, ipa, nitscf, nite, io);
    delete s->vext;
    s->vext = 0;
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
void ResponseCmd::responseVext(bool rpa, bool ipa, int nitscf, int nite, string io)
{
//  if (ui->onpe0())
//  {
//    cout << "ResponseCmd:responseVext:\n"
//         << "RPA = " << rpa << ", IPA = " << ipa << ", io = " << io
//         << ", q = " << q << ", nitscf = " << nitscf << ", nite = " << nite << "\n";
//  }

  s->wf.info(cout, "wavefunction");

  const int nspin = s->wf.nspin();
  BOSampleStepper *stepper = new BOSampleStepper(*s, nitscf, nite);

  if (rpa)
    stepper->set_update_vxc(false);

  if (ipa)
  {
    stepper->set_update_vh(false);
    stepper->set_update_vxc(false);
  }

  // save a copy of initial wave functions
  Wavefunction wf0(s->wf);
  wf0 = s->wf;

  ChargeDensity &cd = stepper->cd();
  MPI_Comm vcomm = cd.vcomm();

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
  const Context* ctxt = s->wf.spincontext();
  int nprow = ctxt->nprow();
  int npcol = ctxt->npcol();
  int myrow = ctxt->myrow();
  int mycol = ctxt->mycol();
  Timer tm_comm_drho, tm_write_drho;

  if (io == "base64_parallel")
  {
#if ! USE_MPI
    cout << "cannot use parallel_io when running in serial" << endl;
#endif
    // MPI parallel writing with base64 encoding
    // all processors on column 1 write collectively
    Base64Transcoder xcdr;
    int lastproc = nprow - 1;
    while ( lastproc >= 0 && ft2.np2_loc(lastproc) == 0 ) lastproc --;
    assert(lastproc >= 0);

    for (int ispin = 0; ispin < nspin; ispin++)
    {
      // following part is a modified version of SlaterDet::write in SlaterDet.C
      // use group of 3 algorithm to adjust the size of drhor on each processor
      tm_comm_drho.start();

      int ndiff;
      assert(drho_r[ispin].size()==n012loc);
      if ( myrow == 0 )
      {
        ndiff = n012loc % 3;
        ctxt->ibcast_send('c',1,1,&ndiff,1);
      }
      else
      {
        ctxt->ibcast_recv('c',1,1,&ndiff,1,0,mycol);
      }
      // send/receive the head/tail of drhor to neighbor
      int nsend_left=0, nsend_right=0, nrecv_left=0, nrecv_right=0;
      if ( myrow % 3 == 0 )
      {
        if ( myrow < lastproc )
          nsend_right = ndiff;
      }
      else if ( myrow % 3 == 1 )
      {
        if ( myrow <= lastproc )
          nrecv_left = ndiff;
        if ( myrow <= lastproc-1 )
          nrecv_right = ndiff;
      }
      else if ( myrow % 3 == 2 )
      {
        if ( myrow <= lastproc && myrow > 0 )
          nsend_left = ndiff;
      }

      double rbuf_left[2], rbuf_right[2], sbuf_left[2], sbuf_right[2];
      int tmpr_size = n012loc;
      if ( nsend_left > 0 )
      {
        for ( int i = 0; i < ndiff; i++ )
          sbuf_left[i] = drho_r[ispin][i];
        ctxt->dsend(ndiff,1,sbuf_left,ndiff,myrow-1,mycol);
        tmpr_size -= ndiff;
      }
      if ( nsend_right > 0 )
      {
        for ( int i = 0; i < ndiff; i++ )
          sbuf_right[i] = drho_r[ispin][n012loc-ndiff+i];
        ctxt->dsend(ndiff,1,sbuf_right,ndiff,myrow+1,mycol);
        tmpr_size -= ndiff;
      }
      if ( nrecv_left > 0 )
      {
        ctxt->drecv(ndiff,1,rbuf_left,ndiff,myrow-1,mycol);
        tmpr_size += ndiff;
      }
      if ( nrecv_right > 0 )
      {
        ctxt->drecv(ndiff,1,rbuf_right,ndiff,myrow+1,mycol);
        tmpr_size += ndiff;
      }

      // at this point, all communications has been finished
      // now construct tmpr to store drhor, tmpr has
      // size divisible by 3 (except last processor)
      if ( myrow < lastproc ) assert(tmpr_size%3 == 0);
      vector<double> tmpr(tmpr_size);
      if ( nrecv_left > 0 || nrecv_right > 0 )
      {
        int index = 0;
        if ( nrecv_left > 0 )
        {
          for ( int i = 0; i < ndiff; i++ )
            tmpr[index++] = rbuf_left[i];
        }
        for ( int i = 0; i < n012loc; i++ )
          tmpr[index++] = drho_r[ispin][i];
        if ( nrecv_right > 0 )
        {
          for ( int i = 0; i < ndiff; i++ )
            tmpr[index++] = rbuf_right[i];
        }
        assert(index==tmpr_size);
      }
      else if ( nsend_left > 0 || nsend_right > 0 )
      {
        int index = 0;
        int istart = (nsend_left > 0) ? ndiff : 0;
        int iend = (nsend_right > 0) ? n012loc - ndiff : n012loc;
        for ( int i = istart; i < iend; i++ )
          tmpr[index++] = drho_r[ispin][i];
        assert(index==tmpr_size);
      }
      else
      {
        for (int i = 0; i < n012loc; i++)
          tmpr[i] = drho_r[ispin][i];
        assert(tmpr_size==n012loc);
      }

      // now transform drhor (as stored in tmpr) as base 64 encoding
      int nbytes = tmpr.size() * sizeof(double);
      int nchars = xcdr.nchars(nbytes);
      char* wbuf = new char[nchars];
      xcdr.encode(nbytes, (byte *) &tmpr[0], wbuf);

      // column 0 write tmpr into file
      string filename;
      if (nspin == 1)
        filename = s->vext->filename() + ".response";
      else
        filename = s->vext->filename() + ".response."
                   + ((ispin == 0) ? "spin0" : "spin1");

      tm_comm_drho.stop();

      if ( mycol == 0 )
      {
        // compute offset
        int offset;
        MPI_Scan(&nchars, &offset, 1, MPI_INT, MPI_SUM, basis.comm());
        offset -= nchars;
        // write file
        tm_write_drho.start();

        MPI_File fh;
        MPI_Info info;
        MPI_Info_create(&info);
        int err;

        err = MPI_File_open(basis.comm(), (char*) filename.c_str(),
                            MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fh);
        if (err != 0)
          cout << myrow << ": error in MPI_File_open: " << err << endl;
        MPI_File_set_size(fh, 0);

        MPI_Status status;
        err = MPI_File_write_at_all(fh,offset,(void*) &wbuf[0],nchars,
                                    MPI_CHAR,&status);
        if ( err != 0 )
          cout << myrow << ": error in MPI_File_write_at_all: err=" << err << endl;

        err = MPI_File_close(&fh);
        if ( err != 0 )
          cout << myrow << ": error in MPI_File_close: " << err << endl;

        tm_write_drho.stop();
      }

      delete [] wbuf;
    } // for ispin
  }
  else if (io == "base64_serial" or io == "cube")
  {
    // serial write with base64 encoding or cube format
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
      tm_comm_drho.start();
      MPI_Gatherv(&drho_r[ispin][0], ft2.np012loc(), MPI_DOUBLE, &drho_r_gathered[0],
                  &rcounts[0], &displs[0], MPI_DOUBLE, 0, vcomm);
      tm_comm_drho.stop();

      if ( myrow == 0 && mycol == 0 )
      {
        string filename;
        if (nspin == 1)
          filename = s->vext->filename() + ".response";
        else
          filename = s->vext->filename() + ".response."
                     + ((ispin == 0) ? "spin0" : "spin1");
        if (io == "base64_serial")
        {
          Base64Transcoder xcdr;
          // transform drhor (stored in drho_r_gathered) to base64 encoding
          int nbytes = drho_r_gathered.size() * sizeof(double);
          int nchars = xcdr.nchars(nbytes);
          char *wbuf = new char[nchars];
          xcdr.encode(nbytes, (byte *) &drho_r_gathered[0], wbuf);

          ofstream os(filename.c_str(), ios::out | ios::binary);
          os.write(wbuf, nchars);
          os.close();
          delete [] wbuf;
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
        } // if io = base64_serial || cube
      } //if ( myrow == 0 && mycol == 0 )
    }
  }
  else
  {
    cout << "unknown io scheme" << endl;
  } // select case io

  double time, tmin, tmax;
  time = tm_comm_drho.real();
  tmin = time;
  tmax = time;
  s->ctxt_.dmin(1, 1, &tmin, 1);
  s->ctxt_.dmax(1, 1, &tmax, 1);
  if (ui->onpe0())
  {
    cout << "  Time to communicate drho "
         << "min: " << tmin << " max: " << tmax << endl;
  }
  time = tm_write_drho.real();
  tmin = time;
  tmax = time;
  s->ctxt_.dmin('C', 1, 1, &tmin, 1);
  s->ctxt_.dmax('C', 1, 1, &tmax, 1);
  if (ui->onpe0())
  {
    cout << "  Time to write drho "
         << "min: " << tmin << " max: " << tmax << endl;
  }

  delete stepper;
}
