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
// SDCellStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "SDCellStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::compute_new_cell(double e0, const valarray<double>& sigma,
  const std::vector<std::vector< double> >& f0)
{
  // compute new cell by adding correction dcell_ij. Steepest descent algorithm.
  // forces f0 are not used
  valarray<double> dcell(9);

  const UnitCell& cell = s_.wf.cell();
  const double cell_mass = s_.ctrl.cell_mass;

  if ( cell_mass <= 0.0 )
  {
    if ( s_.ctxt_.onpe0() )
    {
      cout << " SDCellStepper::compute_new_cell: cell mass is zero\n"
           << "     cannot update cell" << endl;
      return;
    }
  }

  // Compute cell correction dcell_ij from the stress tensor
  // dcell = - omega * sigma * A
  assert(sigma.size()==6);
  // next line: local copy of sigma to circumvent compiler error
  valarray<double> sigma_loc(sigma);
  cell.smatmult3x3(&sigma_loc[0],cell.amat(),&dcell[0]);
  // cell.smatmult3x3(&sigma[0],cell.amat(),&dcell[0]);
  dcell *= -cell.volume();

  // new cell: cellp = cell - dcell * dt^2 / cell_mass
  const double dt = s_.ctrl.dt;
  const double dt2bym = dt*dt/cell_mass;
  D3vector a0p = cell.a(0) - dt2bym * D3vector(dcell[0],dcell[1],dcell[2]);
  D3vector a1p = cell.a(1) - dt2bym * D3vector(dcell[3],dcell[4],dcell[5]);
  D3vector a2p = cell.a(2) - dt2bym * D3vector(dcell[6],dcell[7],dcell[8]);
  cellp = UnitCell(a0p,a1p,a2p);

  // check for cell_lock constraints and modify cellp if needed
  enforce_constraints(cell, cellp);
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
  s_.atoms.set_positions(r);
  s_.atoms.set_cell(cellp);

  // resize wavefunction and basis sets
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
}
