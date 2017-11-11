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
// ConfinementPotential.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CONFINEMENTPOTENTIAL_H
#define CONFINEMENTPOTENTIAL_H

#include <valarray>
class Basis;

class ConfinementPotential
{
  private:

  double ecuts_, facs_, sigmas_;
  const Basis& basis_;
  std::valarray<double> fstress_, dfstress_;

  public:

  double facs(void) const { return facs_; }
  double sigmas(void) const { return sigmas_; }
  double ecuts(void) const { return ecuts_; }

  const std::valarray<double>& fstress(void) const { return fstress_; }
  const std::valarray<double>& dfstress(void) const { return dfstress_; }

  void update(void);

  const Basis& basis() const { return basis_; }

  ConfinementPotential(double ecuts, double facs, double sigmas,
    const Basis& basis);
  ~ConfinementPotential() {}
};
#endif
