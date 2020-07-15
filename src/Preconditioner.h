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
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Wavefunction;
class EnergyFunctional;

#include <vector>
#include <valarray>

class Preconditioner
{
  private:

  EnergyFunctional& ef_;

  // kpg2_[isp_loc][ikp_loc][ig]
  std::vector<std::vector<const double *> > kpg2_;
  // ekin_[isp_loc][ikp_loc][n]
  std::vector<std::vector<std::valarray<double> > > ekin_;

  double ecutprec_;

  public:

  // update values of ekin_
  void update(const Wavefunction& wf);

  double diag(int isp_loc, int ikp_loc, int n, int ig) const;

  Preconditioner(const Wavefunction& wf, EnergyFunctional& ef, double ecutprec);
};
#endif
