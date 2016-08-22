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

int ResponseCmd::action(int argc, char **argv)
{
  // " syntax: response amplitude nitscf [nite]\n\n"
  // " syntax: response -vext vext_file nitscf [nite]\n\n"
  if ( argc < 3 || argc > 5)
  {
    if ( ui->onpe0() )
      cout << " use: response amplitude nitscf [nite]\n"
           << "      response -vext vext_file nitscf [nite]\n"
           << endl;
    return 1;
  }
  if ( s->wf.nst() == 0 )
  {
    if ( ui->onpe0() )
      cout << " ResponseCmd: no states, cannot run" << endl;
    return 1;
  }
  if ( s->wf.ecut() == 0.0 )
  {
    if ( ui->onpe0() )
      cout << " ResponseCmd: ecut = 0.0, cannot run" << endl;
    return 1;
  }

  int nitscf;
  int nite = 0;
  if ( !strcmp(argv[1],"-vext") )
  {
    // response to vext
    if ( ui->onpe0() )
      cout << " ResponseCmd: start computing charge density response under "
           << " external potential from " << argv[2] << endl;
    s->vext->filename = string(argv[2]);
    nitscf = atoi(argv[3]);
    if ( argc == 4 )
      nite = atoi(argv[4]);

    responseVext(nitscf,nite);
    //s->vext->filename.clear();
  }
  else
  {
    // polarizability calculation
    if ( ui->onpe0() )
      cout << " ResponseCmd: start polarizability calculation with "
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

void ResponseCmd::responseVext(int nitscf, int nite)
{

  SampleStepper* stepper = new BOSampleStepper(*s,nitscf,nite);

  stepper->step(0);
  ChargeDensity rho1(s->wf);
  rho1.update_density();
  const vector<vector<double> > &rhor1 = rho1.rhor;

  s->vext->reverse();
  if ( ui->onpe0() )
    cout << " ResponseCmd: external potential is reversed, "
       << " starting another scf iteration" << endl;

  stepper->step(0);
  ChargeDensity rho2(s->wf);
  rho2.update_density();
  const vector<vector<double> > &rhor2 = rho2.rhor;

  // compute drho_r as rhor1 - rhor2
  const int nspin = s->wf.nspin();
  const int np012loc = rho1.vft()->np012loc();
  const double omega = rho1.vbasis()->cell().volume();
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
  const UnitCell& cell = rho1.vbasis()->cell();
  const double ecut = s->vext->ecut();
  const int np0 = rho1.vft()->np0();
  const int np1 = rho1.vft()->np1();
  const int np2 = rho1.vft()->np2();
  Basis basis(rho1.vcomm(),D3vector(0,0,0));
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
        drho_r[ispin][ir] = real( drho_r_tmp[ispin][ir] * omega_inv );
      }
    }

  MPI_Barrier(MPI_COMM_WORLD);
  if ( ui->onpe0() )
    cout << " ResponseCmd: density response has been"
       << " interpolated from grid size \n"
         << " (" << np0 << ", " << np1 << ", " << np2 << ")"
     << " to (" << n0 << ", " << n1 << ", " << n2 << ") \n";

  // process 0 collects the drho_r from all processors on column 0
  const Context& ctxt = *(s->wf.spincontext());
  if ( ctxt.onpe0() )
  {
    for ( int ispin = 0; ispin < nspin; ispin++ )
    drho_r[ispin].resize(ft2.np012(), 0.0);
  }

  if ( ctxt.mycol() == 0 )
  {
  for ( int ispin = 0; ispin < nspin; ispin++ )
  {
    for ( int irow = 0; irow < ctxt.nprow(); irow++ )
    {
    bool iamsending = irow == ctxt.myrow();
    // send drho_r block size
    int size = -1;
    if ( ctxt.onpe0() )
    {
      if ( iamsending )
      {
      }
      else
      {
        ctxt.irecv(1,1,&size,1,irow,0);
        //cout << "I am receiving size from pe " << irow << " size: " << size << endl;
      }
    }// if onpe0
    else
    {
      if ( iamsending )
      {
      size = ft2.np012loc();
      ctxt.isend(1,1,&size,1,0,0);
      //cout << "I am sending size from pe " << ctxt.mype() << " size: " << size << endl;
      }
    }

        // send drho_r block
    if ( ctxt.onpe0() )
    {
      if ( iamsending )
      {
      }
      else
      {
      int istart = ft2.np0() * ft2.np1() * ft2.np2_first(irow);
            ctxt.drecv(size,1,&drho_r[ispin][istart],1,irow,0);
      }
    }
    else
    {
        if ( iamsending )
      {
      ctxt.dsend(size,1,&drho_r[ispin][0],1,0,0);
      }
    }
    } // for irow
    } // for ispin
  } // if mycol = 0

  // process 0 output density difference
  if ( ctxt.onpe0() )
  {
  ofstream os;
    string filename = s->vext->filename + ".response";
  os.open(filename.c_str());
    os << n0 << " " << n1 << " " << n2 << " " << endl;
    for ( int i = 0; i < drho_r[0].size(); i++ )
      if ( nspin == 2 )
        os <<  drho_r[0][i] + drho_r[1][i] << endl;
      else
        os <<  drho_r[0][i] * omega / ft2.np012() << endl;
    cout << " ResponseCmd: charge density response has been written in: "
       << filename << endl;
  }

  delete stepper;
}


