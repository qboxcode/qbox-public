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
// CGCellStepper.C
//
////////////////////////////////////////////////////////////////////////////////

#include "CGCellStepper.h"
#include <valarray>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
CGCellStepper::CGCellStepper(Sample& s) : CellStepper(s), first_step_(true),
  sigma1_(0.1), sigma2_(0.5)
{
  nat_ = atoms_.size();

  xc_.resize(3*nat_+9);
  pc_.resize(3*nat_+9);
  fc_.resize(3*nat_+9);

  rc_.resize(atoms_.nsp());
  rp_.resize(atoms_.nsp());
  for ( int is = 0; is < rc_.size(); is++ )
  {
    rc_[is].resize(3*atoms_.na(is));
    rp_[is].resize(3*atoms_.na(is));
  }

  linmin_.set_sigma1(sigma1_);
}

////////////////////////////////////////////////////////////////////////////////
void CGCellStepper::compute_new_cell(double e0, const valarray<double>& sigma,
  const std::vector<std::vector< double> >& fion)
{
  // compute new cell and ionic positions using the stress tensor sigma
  // and the forces on ions fion
  const UnitCell& cell = s_.wf.cell();

  // total number of dofs: 3* natoms + cell parameters
  valarray<double> f0(3*nat_+9), x0(3*nat_+9), xp(3*nat_+9);

  // copy current positions into x0
  vector<vector<double> > r0;
  atoms_.get_positions(r0);
  double tmp3[3];
  // convert position from r0 to tau (relative) coordinates: tau = A^-1 R
  for ( int i = 0, is = 0; is < r0.size(); is++ )
  {
    for ( int ia = 0; ia < r0[is].size()/3; ia++ )
    {
      cell.vecmult3x3(cell.amat_inv(),&r0[is][3*ia],tmp3);
      x0[i++]=tmp3[0];
      x0[i++]=tmp3[1];
      x0[i++]=tmp3[2];
    }
  }

  // convert forces on positions to forces on tau coordinates, store in f0
  // f = A^-1 * fion
  for ( int i = 0, is = 0; is < fion.size(); is++ )
  {
    for ( int ia = 0; ia < fion[is].size()/3; ia++ )
    {
      cell.vecmult3x3(cell.amat_inv(),&fion[is][3*ia],tmp3);
      f0[i++]=tmp3[0];
      f0[i++]=tmp3[1];
      f0[i++]=tmp3[2];
    }
  }

  // convert descent direction dcell to cell space from the stress tensor
  // dcell = - omega * sigma * A
  valarray<double> dcell(9); // descent direction in cell space
  assert(sigma.size()==6);
  // next line: local copy of sigma to circumvent compiler error
  valarray<double> sigma_loc(sigma);
  cell.smatmult3x3(&sigma_loc[0],cell.amat(),&dcell[0]);
  // dcell *= -cell.volume();
  for ( int i = 0; i < dcell.size(); i++ )
    f0[3*nat_+i] = dcell[i];

  // cout << "dcell = " << endl;
  // cout << dcell[0] << " " << dcell[1] << " " << dcell[2] << endl;
  // cout << dcell[3] << " " << dcell[4] << " " << dcell[5] << endl;
  // cout << dcell[6] << " " << dcell[7] << " " << dcell[8] << endl;

  // the vector f0 now contains the derivatives of the energy in tau+cell space
#ifdef DEBUG
  if ( s_.ctxt_.onpe0() )
  {
    cout << " f0:" << endl;
    for ( int i = 0; i < f0.size(); i++ )
      cout << f0[i] << endl;
  }
#endif

  double fp0; // derivative f'(x) at x=xp
  bool wolfe1, wolfe2; // First and second Wolfe conditions in line search

  // CG algorithm

  if ( !first_step_ )
  {
    wolfe1 = e0 < ec_ + fpc_ * sigma1_ * alpha_;
    // compute fp0: projection of forces on direction pc_
    fp0 = 0.0;
    for ( int i = 0; i < f0.size(); i++ )
      fp0 -= f0[i] * pc_[i];

    wolfe2 = fabs(fp0) < sigma2_ * fabs(fpc_);
    if ( s_.ctxt_.onpe0() )
    {
      cout << "  CGCellStepper: fpc = " << fpc_ << endl;
      cout << "  CGCellStepper: fp0 = " << fp0 << endl;
      cout << "  CGCellStepper: ec = " << ec_ << " e0 = " << e0 <<  endl;
      cout << "  CGCellStepper: ec_ + fpc_ * sigma1_ * alpha_ ="
           << ec_ + fpc_ * sigma1_ * alpha_ << endl;
      cout << "  CGCellStepper: wolfe1/wolfe2 = "
           << wolfe1 << "/" << wolfe2 << endl;
    }
  }

  if ( first_step_ || (wolfe1 && wolfe2) )
  {
    // define new descent direction pc_
    if ( first_step_ )
    {
      pc_ = f0;
    }
    else
    {
      // Polak-Ribiere definition
      double num = 0.0, den = 0.0;
      for ( int i = 0; i < f0.size(); i++ )
      {
        const double fctmp = fc_[i];
        const double f0tmp = f0[i];
        num += f0tmp * ( f0tmp - fctmp );
        den += fctmp * fctmp;
      }
      double beta = den > 0.0 ? num/den : 0.0;
      beta = max(beta,0.0);
      if ( s_.ctxt_.onpe0() )
        cout << "  CGCellStepper: beta = " << beta << endl;
      for ( int i = 0; i < f0.size(); i++ )
        pc_[i] = beta * pc_[i] + f0[i];
    }

    fc_ = f0;
    // fpc = d_e / d_alpha in direction pc
    fpc_ = 0.0;
    for ( int i = 0; i < f0.size(); i++ )
      fpc_ -= fc_[i] * pc_[i];
    ec_ = e0;
    xc_ = x0;
    fp0 = fpc_;
    // reset line minimizer
    linmin_.reset();
  }

#ifdef DEBUG
  if ( s_.ctxt_.onpe0() )
  {
    cout << " pc:" << endl;
    for ( int i = 0; i < pc_.size(); i++ )
      cout << pc_[i] << endl;
  }

  // find largest component of pc_
  double pcmax = 0.0;
  for ( int i = 0; i < pc_.size(); i++ )
    pcmax = max(fabs(pc_[i]),pcmax);
  if ( s_.ctxt_.onpe0() )
  {
    cout << "CGCellStepper: largest component of pc_: "
         << pcmax << endl;
  }
#endif

  alpha_ = linmin_.newalpha(alpha_,e0,fp0);

  if ( s_.ctxt_.onpe0() )
    cout << "  CGCellStepper: alpha = " << alpha_ << endl;

  // update
  xp = xc_ + alpha_ * pc_;

  // save the change of cell parameters into dcell
  for ( int i = 0; i < 9; i++ )
    dcell[i] = xp[3*nat_+i]-xc_[3*nat_+i];

  // cellp = cell + dcell
  D3vector a0p = cell.a(0) + D3vector(dcell[0],dcell[1],dcell[2]);
  D3vector a1p = cell.a(1) + D3vector(dcell[3],dcell[4],dcell[5]);
  D3vector a2p = cell.a(2) + D3vector(dcell[6],dcell[7],dcell[8]);
  // cout << "dcell = " << endl;
  // cout << dcell[0] << " " << dcell[1] << " " << dcell[2] << endl;
  // cout << dcell[3] << " " << dcell[4] << " " << dcell[5] << endl;
  // cout << dcell[6] << " " << dcell[7] << " " << dcell[8] << endl;
  cellp = UnitCell(a0p,a1p,a2p);

  // check for cell_lock constraints and modify cellp if needed
  enforce_constraints(cell, cellp);

  // compute atomic positions rc_[is][3*ia+j] from xc_
  // compute new atomic positions rp_[is][3*ia+j] from xp and cellp
  for ( int is = 0, i = 0; is < fion.size(); is++ )
  {
    for ( int ia = 0; ia < fion[is].size()/3; ia++ )
    {
      tmp3[0] = xc_[i+0];
      tmp3[1] = xc_[i+1];
      tmp3[2] = xc_[i+2];
      cell.vecmult3x3(cell.amat(),tmp3,&rc_[is][3*ia]);
      tmp3[0] = xp[i+0];
      tmp3[1] = xp[i+1];
      tmp3[2] = xp[i+2];
      cell.vecmult3x3(cellp.amat(),tmp3,&rp_[is][3*ia]);
      i+=3;
    }
  }

  // enforce constraints on atomic positions
  s_.constraints.enforce_r(rc_,rp_);

  first_step_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void CGCellStepper::update_cell(void)
{
  const UnitCell& cell = s_.wf.cell();

  s_.atoms.set_positions(rp_);
  s_.atoms.set_cell(cellp);

  // resize wavefunction and basis sets
  s_.wf.resize(cellp,s_.wf.refcell(),s_.wf.ecut());
  if ( s_.wfv != 0 )
    s_.wfv->resize(cellp,s_.wf.refcell(),s_.wf.ecut());
}
