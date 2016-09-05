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
// XCOperator.C
//
////////////////////////////////////////////////////////////////////////////////
#include "XCOperator.h"
#include "ChargeDensity.h"
#include "XCPotential.h"
#include "ExchangeOperator.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
XCOperator::XCOperator(Sample& s, const ChargeDensity& cd) :cd_(cd)
{
  // set initial values
  xcp_ = 0;
  xop_ = 0;
  exc_ = 0.0 ;
  dxc_ = 0.0 ;

  sigma_exc_.resize(6);

  string functional_name = s.ctrl.xc;

  // check the name of the functional
  if ( ( functional_name ==  "LDA" ) ||
       ( functional_name ==  "VWN" ) ||
       ( functional_name ==  "PBE" ) ||
       ( functional_name == "BLYP" ) )
  {
    // create only an xc potential
    xcp_ = new XCPotential(cd, functional_name, s.ctrl);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = false;
    HFmixCoeff_ = 0.0;
  }
  else if ( functional_name == "HF" )
  {
    // create exchange operator with mixing coeff=1
    xop_ = new ExchangeOperator(s, 1.0);
    hasPotential_ = false;
    hasGGA_ = false;
    hasHF_ = true;
    HFmixCoeff_ = 1.0;
  }
  else if ( functional_name == "PBE0" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s.ctrl);

    // create the exchange operator with mixing coeff=0.25
    xop_ = new ExchangeOperator(s, s.ctrl.alpha_PBE0);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
    HFmixCoeff_ = s.ctrl.alpha_PBE0;;
  }
  else if ( functional_name == "B3LYP" )
  {
    // create an exchange potential
    xcp_ = new XCPotential(cd, functional_name, s.ctrl);

    // create the exchange operator with mixing coeff=0.20
    xop_ = new ExchangeOperator(s, 0.20);
    hasPotential_ = true;
    hasGGA_ = xcp_->isGGA();
    hasHF_ = true;
    HFmixCoeff_ = 0.20;
  }
  else
  {
    throw XCOperatorException(
      "unknown functional name during exchange operator construction");
  }
}

XCOperator::~XCOperator()
{
  delete xcp_;
  delete xop_;
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::update(std::vector<std::vector<double> >& vr,
                        bool compute_stress, bool freeze_vxc)
{
  // update xc potential and self-energy
  // used whenever the charge density and/or wave functions have changed
  // compute vxc potential and energy
  if ( hasPotential_ )
  {
    // update LDA/GGA xc potential
    xcp_->update(vr, freeze_vxc);

    // LDA/GGA exchange energy
    exc_ = xcp_->exc();
    dxc_ = xcp_->dxc();

    if ( compute_stress )
      xcp_->compute_stress(sigma_exc_);
  }
  else
  {
    exc_ = 0.0;
    dxc_ = 0.0;
    sigma_exc_ = 0.0;
  }

  if ( hasHF() )
  {
    double ex_hf = xop_->update_operator(compute_stress);
    exc_ += ex_hf;
    dxc_ -= ex_hf;
    if ( compute_stress )
      xop_->add_stress(sigma_exc_);
  }
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::apply_self_energy(Wavefunction &dwf)
{
  if ( hasHF() )
    xop_->apply_operator(dwf);
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::compute_stress(std::valarray<double>& sigma)
{
  sigma = sigma_exc_;
}

////////////////////////////////////////////////////////////////////////////////
void XCOperator::cell_moved(void)
{
  if ( hasHF() )
    xop_->cell_moved();
}
