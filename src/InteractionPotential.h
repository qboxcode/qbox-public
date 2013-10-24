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
// InteractionPotential.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Implements InteractionPotential class that evaluates the potential for given
// norm of G vectors and the derivative w.r.t. this argument
//
// Author: Martin Schlipf (2013)
// Contact: martin.schlipf@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef INTERACTIONPOTENTIAL_H
#define INTERACTIONPOTENTIAL_H

#include <cassert>

class InteractionPotential
{
  public:

  // default constructor = Coulomb potential
  InteractionPotential() :
    coulomb_(true)
  {
  }

  // constructor - define function and derivative
  InteractionPotential(double(*V)(const double&), double(*dV)(const double&),
    double(*div_scale)(const double&)) :
    V_(V), dV_(dV), div_scale_(div_scale), coulomb_(false)
  {
  }

  // is the interaction potential a coulomb potential?
  inline bool coulomb() const
  {
    return coulomb_;
  }

  // evaluate the interaction potential for given norm of G vector
  inline double operator()(const double G2) const
  {
    // the current implementation expects that coulomb potential is treated externaly
    assert(not coulomb_);
    return V_(G2);
  }

  // evaluate the derivative of the interaction potential w.r.t. G^2
  inline double derivative(const double G2) const
  {
    // the current implementation expects that coulomb potential is treated externaly
    assert(not coulomb_);
    return dV_(G2);
  }

  inline double divergence_scaling(const double rcut) const
  {
    // the current implementation expects that coulomb potential is treated externaly
    assert(not coulomb_);
    return div_scale_(rcut);
  }

  private:

  const bool coulomb_;
  double (*V_)(const double&);
  double (*dV_)(const double&);
  double (*div_scale_)(const double&);

};

#endif
