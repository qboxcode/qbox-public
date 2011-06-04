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
// CellStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "CellStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void CellStepper::enforce_constraints(const UnitCell& cell, UnitCell& cellp)
{
  // modify cellp so as to satisfy the constraints expressed in cell_lock
  string cell_lock = s_.ctrl.cell_lock;

  valarray<double> dcell(9); // cellp - cell
  for ( int i = 0; i < 9; i++ )
    dcell[i] = cellp.amat(i) - cell.amat(i);

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

  // dcell now represents a change compatible with constraints
  // cellp = cell + dcell
  D3vector a0p = cell.a(0) + D3vector(dcell[0],dcell[1],dcell[2]);
  D3vector a1p = cell.a(1) + D3vector(dcell[3],dcell[4],dcell[5]);
  D3vector a2p = cell.a(2) + D3vector(dcell[6],dcell[7],dcell[8]);
  cellp = UnitCell(a0p,a1p,a2p);
}
