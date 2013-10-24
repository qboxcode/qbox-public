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

class InteractionPotential
{
  public:

  // constructor - define function and derivative
  InteractionPotential(double(*V)(const double&), double(*dV)(const double&)) :
    V_(V), dV_(dV)
  {
  }

  // evaluate the interaction potential for given norm of G vector
  inline double operator()(const double G2) const
  {
    return V_(G2);
  }

  // evaluate the derivative of the interaction potential w.r.t. G^2
  inline double derivative(const double G2) const
  {
    return dV_(G2);
  }

  private:

  double (*V_)(const double&);
  double (*dV_)(const double&);

};

#endif
