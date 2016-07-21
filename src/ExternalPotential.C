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
// ExternalPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include "ExternalPotential.h"
#include "FourierTransform.h"
#include "Basis.h"
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

void ExternalPotential::update(const ChargeDensity& cd)
{
  // read vext on nrow=0 tasks
  vector<double> vext_read;
  FourierTransform *vft = cd.vft();

  const int myrow = cd.context().myrow();
  if ( myrow == 0 )
  {
    // read the external potential from file s_.ctrl.vext
    if ( cd.context().onpe0() )
      cout << "ExternalPotential::update: vext file: " << s_.ctrl.vext << endl;
    ifstream vfile(s_.ctrl.vext.c_str());
    vfile >> n_[0] >> n_[1] >> n_[2];
    if ( cd.context().onpe0() )
      cout << "ExternalPotential::update: grid size "
           << n_[0] << " " << n_[1] << " " << n_[2] << endl;
    int n012 = n_[0] * n_[1] * n_[2];
    vext_read.resize(n012);
    for ( int i = 0; i < n012; i++ )
      vfile >> vext_read[i];
  }
  // broadcast sizes to lower rows
  MPI_Bcast(&n_[0],3,MPI_INT,0,cd.vcomm());
  int n012 = n_[0] * n_[1] * n_[2];
  vext_read.resize(n012);
  MPI_Bcast(&vext_read[0],n012,MPI_DOUBLE,0,cd.vcomm());
  // vext_read now contains vext data on all tasks

  // todo: replace Bcast with Scatter

  // create a Basis compatible with the vext grid read from file
  // determine largest ecut compatible with grid and current unit cell
  const Basis* vbasis = cd.vbasis();
  const UnitCell& cell = vbasis->cell();
  const D3vector b0 = cell.b(0);
  const D3vector b1 = cell.b(1);
  const D3vector b2 = cell.b(2);

  int n0max = (n_[0]-2)/2;
  int n1max = (n_[1]-2)/2;
  int n2max = (n_[2]-2)/2;
  double ecut0 = 0.5 * norm2(b0) * n0max*n0max;
  double ecut1 = 0.5 * norm2(b1) * n1max*n1max;
  double ecut2 = 0.5 * norm2(b2) * n2max*n2max;

  double ecut = min(min(ecut0,ecut1),ecut2);

  if ( cd.context().onpe0() )
  {
    cout << "ExternalPotential::update: ecut0: " << ecut0 << endl;
    cout << "ExternalPotential::update: ecut1: " << ecut1 << endl;
    cout << "ExternalPotential::update: ecut2: " << ecut2 << endl;
    cout << "ExternalPotential::update: ecut:  " << ecut << endl;
  }

  Basis basis(cd.vcomm(),D3vector(0,0,0));
  basis.resize(cell,cell,ecut);
  if ( cd.context().onpe0() )
  {
    cout << "ExternalPotential::update: np0: " << basis.np(0) << endl;
    cout << "ExternalPotential::update: np1: " << basis.np(1) << endl;
    cout << "ExternalPotential::update: np2: " << basis.np(2) << endl;
  }

  // interpolate on grid compatible with the charge density cd
  FourierTransform ft(basis,n_[0],n_[1],n_[2]);
  vector<complex<double> > vext_read_g(basis.localsize());
  vector<complex<double> > tmp_r(ft.np012loc());
  // index of local vext slice in global vext array
  int istart = ft.np2_first() * n_[0] * n_[1];
  for ( int i = 0; i < tmp_r.size(); i++ )
    tmp_r[i] = vext_read[istart+i];

  // compute Fourier coefficients
  ft.forward(&tmp_r[0],&vext_read_g[0]);
  // vext_read_g now contains the Fourier coefficients of vext
  // interpolate to vft grid
  FourierTransform ft2(basis,vft->np0(),vft->np1(),vft->np2());
  tmp_r.resize(vft->np012loc());
  vext_r_.resize(tmp_r.size());
  ft2.backward(&vext_read_g[0],&tmp_r[0]);
  for ( int i = 0; i < vext_r_.size(); i++ )
    vext_r_[i] = real(tmp_r[i]);
}
