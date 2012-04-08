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
// XCOperator.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef XCOPERATOR_H
#define XCOPERATOR_H

#include "Sample.h"
#include "ChargeDensity.h"

class XCPotential;
class ExchangeOperator;
class XCOperator
{
  private:

  XCPotential* xcp_;
  ExchangeOperator* xop_;

  const ChargeDensity& cd_;
  double HFmixCoeff_ ;
  double exc_; // local XC energy
  double exhf_; // Hartree-Fock exchange energy

  bool hasPotential_;
  bool hasGGA_;
  bool hasHF_;

  public:

  // constructor
  XCOperator( Sample& s, const ChargeDensity& cd);

  // destructor
  ~XCOperator();

  // return pointer to the exchange potential class
  const XCPotential* xcp() { return xcp_; }

  // return pointer to the ExchangeOperator class
  ExchangeOperator* xop() { return xop_; }

  bool hasGGA(void) { return hasGGA_; };
  bool hasHF(void) { return hasHF_; };

  void update_v(std::vector<std::vector<double> >& vr);
  double update_sigma(void);
  void apply_sigma(Wavefunction &dwf);
  void compute_stress(std::valarray<double>& sigma_exc);
  double exc(void) { return exc_; };
  double exhf(void) { return exhf_; };
};

class XCOperatorException
{
  public:
  std::string msg;
  XCOperatorException(std::string s) : msg(s) {}
};

#endif
