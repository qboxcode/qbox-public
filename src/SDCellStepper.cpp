////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// SDCellStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "SDCellStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SDCellStepper::SDCellStepper(Sample& s) : CellStepper(s)
{
  // strain tensors u, up
  u_.resize(9);
  up_.resize(9);
  for ( int i = 0; i < u_.size(); i++ )
  {
    u_[i] = 0.0;
    up_[i] = 0.0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::compute_new_cell(double e0, const valarray<double>& sigma,
  const std::vector<std::vector< double> >& fion)
{
  // compute new cell and ionic positions using the stress tensor sigma
  const UnitCell cell0 = s_.wf.cell();
  const double cell_mass = s_.ctrl.cell_mass;

  if ( cell_mass <= 0.0 )
  {
    cellp = cell0;
    if ( MPIdata::onpe0() )
    {
      cout << " SDCellStepper::compute_new_cell: cell mass is zero\n"
           << "     cannot update cell" << endl;
      return;
    }
  }

  assert(sigma.size()==6);
  const double dt = s_.ctrl.dt;
  const double dt2bym = dt*dt/cell_mass;

  double g[9];
  // compute descent direction in strain space
  const double omega = cell0.volume();
  // gradient g = - sigma * volume
  g[0] = -sigma[0] * omega;
  g[1] = -sigma[3] * omega;
  g[2] = -sigma[5] * omega;

  g[3] = -sigma[3] * omega;
  g[4] = -sigma[1] * omega;
  g[5] = -sigma[4] * omega;

  g[6] = -sigma[5] * omega;
  g[7] = -sigma[4] * omega;
  g[8] = -sigma[2] * omega;

  // the vector g now contains the gradient of the energy in strain space

  // SD algorithm
  for ( int i = 0; i < 9; i++ )
    up_[i] = - dt2bym * g[i];

  // check for cell_lock constraints and modify up_ if needed
  enforce_constraints(&up_[0]);

  // compute cellp = ( I + Up ) * A0
  // symmetric matrix I+U stored in iumat: xx, yy, zz, xy, yz, xz
  double iupmat[6];
  iupmat[0] = 1.0 + up_[0];
  iupmat[1] = 1.0 + up_[4];
  iupmat[2] = 1.0 + up_[8];
  iupmat[3] = 0.5 * ( up_[1] + up_[3] );
  iupmat[4] = 0.5 * ( up_[5] + up_[7] );
  iupmat[5] = 0.5 * ( up_[2] + up_[6] );

  const double *a0mat = cell0.amat();
  double apmat[9];
  cell0.smatmult3x3(&iupmat[0],a0mat,apmat);

  D3vector a0p(apmat[0],apmat[1],apmat[2]);
  D3vector a1p(apmat[3],apmat[4],apmat[5]);
  D3vector a2p(apmat[6],apmat[7],apmat[8]);
  cellp.set(a0p,a1p,a2p);
}

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::update_cell(void)
{
  const UnitCell& cell = s_.wf.cell();

  // rescale atomic positions in AtomSet

  // r_new = A_new A_old^-1 r_old
  vector<vector<double> > r;
  s_.atoms.get_positions(r);

  double tau[3];
  for ( int is = 0; is < r.size(); is++ )
  {
    // transform r to tau: multiply by A^-1
    const int nais = r[is].size()/3;
    for ( int ia = 0; ia < nais; ia++ )
    {
      // multiply r[is][ia] by A_old^-1, result in tau
      cell.vecmult3x3(cell.amat_inv(),&r[is][3*ia],&tau[0]);
      // multiply tau by A_new, result in r[is][3*ia]
      cellp.vecmult3x3(cellp.amat(),&tau[0],&r[is][3*ia]);
    }
  }
  s_.atoms.sync_positions(r);
  s_.atoms.set_positions(r);
  s_.atoms.sync_cell(cellp);
  s_.atoms.set_cell(cellp);

  // resize wavefunction and basis sets
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
}
