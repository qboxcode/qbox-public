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
// $Id: SDCellStepper.C,v 1.11 2008-09-26 21:05:18 fgygi Exp $

#include "SDCellStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SDCellStepper::compute_new_cell(const valarray<double>& sigma)
{
  // compute new cell by adding correction dcell_ij
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

  string cell_lock = s_.ctrl.cell_lock;
  if ( cell_lock != "OFF" )
  {
    // constraints on the cell derivatives
    if ( cell_lock.find("A") != string::npos )
    {
      // vector A is locked
      dcell[0] = dcell[1] = dcell[2] = 0.0;
    }
    if ( cell_lock.find("B") != string::npos )
    {
      // vector B is locked
      dcell[3] = dcell[4] = dcell[5] = 0.0;
    }
    if ( cell_lock.find("C") != string::npos )
    {
      // vector C is locked
      dcell[6] = dcell[7] = dcell[8] = 0.0;
    }

    // Check if cell shape should be preserved (if "S" is present in cell_lock)
    // The only changes allowed are renormalizations of a,b,c
    if ( cell_lock.find("S") != string::npos )
    {
      // projection of d in the direction of e
      D3vector d,e;

      d = D3vector(dcell[0],dcell[1],dcell[2]);
      e = cell.a(0) / length(cell.a(0));
      d = (d * e) * e;
      dcell[0] = d.x; dcell[1] = d.y; dcell[2] = d.z;

      d = D3vector(dcell[3],dcell[4],dcell[5]);
      e = cell.a(1) / length(cell.a(1));
      d = (d * e) * e;
      dcell[3] = d.x; dcell[4] = d.y; dcell[5] = d.z;

      d = D3vector(dcell[6],dcell[7],dcell[8]);
      e = cell.a(2) / length(cell.a(2));
      d = (d * e) * e;
      dcell[6] = d.x; dcell[7] = d.y; dcell[8] = d.z;
    }

    if ( cell_lock == "R" )
    {
      // preserve aspect ratio
      // dcell must be proportional to A
      // All vectors are rescaled by the same constant
      // rescale cell by coefficient alpha, i.e.
      // dcell = alpha * A
      // where alpha * A is the projection of dcell in the direction of A
      // alpha = tr (A^T * dcell) / || A ||^2  (matrix scalar product)
      const double *a = cell.amat();
      const double num = a[0]*dcell[0] + a[1]*dcell[1] + a[2]*dcell[2] +
                         a[3]*dcell[3] + a[4]*dcell[4] + a[5]*dcell[5] +
                         a[6]*dcell[6] + a[7]*dcell[7] + a[8]*dcell[8];
      const double denom = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] +
                           a[3]*a[3] + a[4]*a[4] + a[5]*a[5] +
                           a[6]*a[6] + a[7]*a[7] + a[8]*a[8];
      const double alpha = num / denom;
      dcell = valarray<double>(a,9);
      dcell *= alpha;
    }
  }

  const double dt = s_.ctrl.dt;
  const double dt2bym = dt*dt/cell_mass;

  // cellp = cell - dcell * dt^2 / cell_mass
  D3vector a0p = cell.a(0) - dt2bym * D3vector(dcell[0],dcell[1],dcell[2]);
  D3vector a1p = cell.a(1) - dt2bym * D3vector(dcell[3],dcell[4],dcell[5]);
  D3vector a2p = cell.a(2) - dt2bym * D3vector(dcell[6],dcell[7],dcell[8]);

  cellp = UnitCell(a0p,a1p,a2p);
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

  //cout << " SDCellStepper::update_cell" << endl;
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
}
