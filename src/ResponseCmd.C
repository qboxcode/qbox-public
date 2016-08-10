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

#include "ResponseCmd.h"
#include<iostream>
#include<sstream>
using namespace std;
#include "BOSampleStepper.h"

#include<ctime>
#include<cassert>

int ResponseCmd::action(int argc, char **argv)
{
  // " syntax: response amplitude nitscf [nite]\n\n"
  // " syntax: response -vext vext_file nitscf [nite]\n\n"
  if ( argc < 2 || argc > 5)
  {
    if ( ui->onpe0() )
      cout << " use: response amplitude nitscf [nite]" << endl
           << "      response [-vext vext_file] nitscf [nite]" << endl;
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

  double amplitude = 0.0;
  bool vext = false;
  string vext_filename;
  int nitscf = 0;
  int nite = 0;
  if ( !strcmp(argv[1],"-vext") )
  {
    assert(argc==4 || argc==5);
    vext = true;
    vext_filename = atoi(argv[2]);
    nitscf = atoi(argv[3]);
    if ( argc == 5 )
      nite = atoi(argv[4]);
    if ( ui->onpe0() )
      cout << " ResponseCmd: -vext option not implemented" << endl;
    return 1;
  }
  else
  {
    assert(argc==3 || argc==4);
    amplitude = atof(argv[1]);
    if ( amplitude == 0.0 )
    {
      if ( ui->onpe0() )
        cout << " ResponseCmd: amplitude is 0.0" << endl;
      return 1;
    }
    nitscf = atoi(argv[2]);
    if ( argc == 4 )
      nite = atoi(argv[3]);
  }

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

  D3vector ddx = 0.5 * (dipole_p[0]-dipole_m[0])/amplitude;
  D3vector ddy = 0.5 * (dipole_p[1]-dipole_m[1])/amplitude;
  D3vector ddz = 0.5 * (dipole_p[2]-dipole_m[2])/amplitude;

  const UnitCell& cell = s->wf.cell();
  cell.fold_in_ws(ddx);
  cell.fold_in_ws(ddy);
  cell.fold_in_ws(ddz);

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

  return 0;
}
