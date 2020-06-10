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
// CellStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "CellStepper.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void CellStepper::enforce_constraints(double *u)
{
  // modify the strain u so as to satisfy the constraints expressed in cell_lock
  string cell_lock = s_.ctrl.cell_lock;

  if ( cell_lock != "OFF" )
  {
    // constraints on the cell derivatives
    if ( cell_lock.find("A") != string::npos )
    {
      // vector A is locked
      u[0] = 0.0;
    }
    if ( cell_lock.find("B") != string::npos )
    {
      // vector B is locked
      u[4] = 0.0;
    }
    if ( cell_lock.find("C") != string::npos )
    {
      // vector C is locked
      u[8] = 0.0;
    }

    // Check if cell shape should be preserved (if "S" is present in cell_lock)
    // The only changes allowed are renormalizations of a,b,c
    if ( cell_lock.find("S") != string::npos )
    {
      // set all off-diagonal elements of u to zero
      u[1] = u[2] = u[3] = u[5] = u[6] = u[7] = 0.0;
    }

    if ( cell_lock == "R" )
    {
      // set all off-diagonal elements of u to zero
      u[1] = u[2] = u[3] = u[5] = u[6] = u[7] = 0.0;
      // set all diagonal elements to 1/3 trace u
      const double tr3rd = ( u[0] + u[4] + u[8] ) / 3.0;
      u[0] = u[4] = u[8] = tr3rd;
    }
  }

  // u is a strain tensor respecting constraints
}
