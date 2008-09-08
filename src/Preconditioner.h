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
// $Id: Preconditioner.h,v 1.5 2008-09-08 15:56:19 fgygi Exp $

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Sample;
class EnergyFunctional;

#include <vector>
#include <valarray>

class Preconditioner
{
  private:

  const Sample& s_;
  const EnergyFunctional& ef_;
  std::vector<std::vector<std::valarray<double> > > diag_; // diag_[ispin][ikp][ig]

  public:

  void update(void);

  const std::valarray<double>& diag(int ispin, int ikp) const
  { return diag_[ispin][ikp]; }

  Preconditioner(const Sample& s, const EnergyFunctional& ef);
  //~Preconditioner();
};
#endif
