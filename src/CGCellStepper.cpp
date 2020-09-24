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
// CGCellStepper.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "CGCellStepper.h"
#include "CGOptimizer.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CGCellStepper::CGCellStepper(Sample& s) : CellStepper(s),
  cgopt_(CGOptimizer(3*s.atoms.size()+9)), cell0(s_.atoms.cell())
{
  nat_ = atoms_.size();
  cgopt_.set_alpha_start(0.002);
  cgopt_.set_alpha_max(0.5);
  cgopt_.set_beta_max(10.0);
#ifdef DEBUG
  if ( MPIdata::onpe0() )
    cgopt_.set_debug_print();
#endif

  rp_.resize(s.atoms.nsp());
  for ( int is = 0; is < rp_.size(); is++ )
    rp_[is].resize(3*atoms_.na(is));

  // store full strain tensor u for consistency of dot products
  u_.resize(9);
  up_.resize(9);
  for ( int i = 0; i < u_.size(); i++ )
  {
    u_[i] = 0.0;
    up_[i] = 0.0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void CGCellStepper::compute_new_cell(double e, const valarray<double>& sigma,
  const std::vector<std::vector< double> >& fion)
{
  // compute new cell and ionic positions using the stress tensor sigma
  // and the forces on ions fion
  const UnitCell cell = s_.wf.cell();

  // total number of dofs: 3* natoms + cell parameters
  valarray<double> x(3*nat_+9), xp(3*nat_+9), g(3*nat_+9);

  // copy current positions into x
  vector<vector<double> > r0, gvec;
  atoms_.get_positions(r0);
  double tmp3[3];
  // convert position from r0 to tau (relative) coordinates: tau = A^-1 R
  for ( int is = 0, i = 0; is < r0.size(); is++ )
  {
    for ( int ia = 0; ia < r0[is].size()/3; ia++ )
    {
      cell.vecmult3x3(cell.amat_inv(),&r0[is][3*ia],tmp3);
      x[i++]=tmp3[0];
      x[i++]=tmp3[1];
      x[i++]=tmp3[2];
    }
  }

  // copy current strain tensor into x
  for ( int i = 0; i < 9; i++ )
    x[3*nat_+i] = u_[i];

  // convert forces on positions to forces on tau coordinates, store -f in gvec
  // f = A^-1 * fion
  gvec.resize(r0.size());
  for ( int is = 0; is < r0.size(); is++ )
  {
    gvec[is].resize(r0[is].size());
    for ( int ia = 0; ia < r0[is].size()/3; ia++ )
    {
      cell.vecmult3x3(cell.amat_inv(),&fion[is][3*ia],tmp3);
      gvec[is][3*ia+0]=-tmp3[0];
      gvec[is][3*ia+1]=-tmp3[1];
      gvec[is][3*ia+2]=-tmp3[2];
    }
  }

  // project the gradient gvec in a direction compatible with constraints
  s_.constraints.enforce_v(r0,gvec);

  for ( int j = 0, is = 0; is < r0.size(); is++ )
    for ( int i = 0; i < r0[is].size(); i++ )
      g[j++] = gvec[is][i];

  // compute descent direction in strain space
  // gradient g = - sigma * volume
  g[3*nat_+0] = -sigma[0] * cell.volume();
  g[3*nat_+1] = -sigma[3] * cell.volume();
  g[3*nat_+2] = -sigma[5] * cell.volume();

  g[3*nat_+3] = -sigma[3] * cell.volume();
  g[3*nat_+4] = -sigma[1] * cell.volume();
  g[3*nat_+5] = -sigma[4] * cell.volume();

  g[3*nat_+6] = -sigma[5] * cell.volume();
  g[3*nat_+7] = -sigma[4] * cell.volume();
  g[3*nat_+8] = -sigma[2] * cell.volume();

  // the vector g now contains the gradient of the energy in tau+strain space

  // enforce constraints on the unit cell
  // Next line only affects the cell components of g
  enforce_constraints(&g[3*nat_]);

  // the vector g now contains a gradient of the energy in tau+strain space
  // compatible with atoms and cell constraints

  // CG algorithm

  cgopt_.compute_xp(x,e,g,xp);

  if ( MPIdata::onpe0() )
  {
    cout << "  CGCellStepper: alpha = " << cgopt_.alpha() << endl;
  }

  for ( int i = 0; i < 9; i++ )
    up_[i] = xp[3*nat_+i];

  // enforce cell_lock constraints
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

  // compute new atomic positions rp_[is][3*ia+j] from xp and cellp
  for ( int is = 0, i = 0; is < fion.size(); is++ )
  {
    for ( int ia = 0; ia < fion[is].size()/3; ia++ )
    {
      tmp3[0] = xp[i+0];
      tmp3[1] = xp[i+1];
      tmp3[2] = xp[i+2];
      cellp.vecmult3x3(cellp.amat(),tmp3,&rp_[is][3*ia]);
      i+=3;
    }
  }

  // enforce constraints on atomic positions
  s_.constraints.enforce_r(r0,rp_);
}

////////////////////////////////////////////////////////////////////////////////
void CGCellStepper::update_cell(void)
{
  s_.atoms.set_positions(rp_);
  s_.atoms.set_cell(cellp);
  u_ = up_;

  // resize wavefunction and basis sets
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
}
