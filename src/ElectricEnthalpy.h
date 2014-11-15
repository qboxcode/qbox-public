////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// ElectricEnthalpy.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ELECTRICENTHALPY_H
#define ELECTRICENTHALPY_H

#include <vector>
#include <complex>
#include <cassert>
#include <string>
#include "Matrix.h"
#include "D3vector.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Context.h"
#include "Sample.h"
#include "Timer.h"
#include "D3tensor.h"
#include "Basis.h"

class Sample;
class MLWFTransform;
class ComplexMatrix;

class ElectricEnthalpy
{
  private:

  Sample& s_;
  Wavefunction& wf_;
  Wavefunction dwf_;
  SlaterDet& sd_;
  const Context& ctxt_;
  const Context& vctxt_;

  bool onpe0_;
  Basis* vbasis_;

  std::string pol_type_;
  int niter_;
  bool compute_quadrupole_;

  // electric field
  D3vector e_field_;

  // phase
  double gamma_[3];

  Wavefunction* rwf_[3];

  // MLWFtransform is used to compute S matrix
  MLWFTransform* mlwft_;

  // s matrices
  ComplexMatrix* smat_[3];

  // total, ionic and electronic part of macroscopic polarization
  D3vector polarization_, polarization_ion_, polarization_elec_;
  D3vector polarization_elec_correction_;
  D3vector polarization_elec_correction_real_;

  // total enthalpy
  double energy_;

  std::vector <D3vector> mlwfc_;
  std::vector <double> mlwfs_;
  std::vector <D3vector> correction_;
  std::vector <D3vector> correction_real_;
  std::vector <D3tensor> quad_;

  void correction(void);
  void correction_real(void);

  void compute_polarization_ion(void);
  void compute_polarization_elec(void);
  void compute_polarization(void);

  public:

  mutable TimerMap tmap;

  bool compute_quadrupole(void) {return compute_quadrupole_;}

  D3vector e_field(){return e_field_;};
  D3vector& polarization(){return polarization_;};
  D3vector& polarization_ion(){return polarization_ion_;};
  D3vector& polarization_elec(){return polarization_elec_;};
  double energy(Wavefunction& dwf, bool compute_hpsi);

  std::vector<D3vector>& mlwfc(void){return mlwfc_;}
  std::vector<double>& mlwfs(void){return mlwfs_;}
  std::vector<D3vector>& cor_reci(void){return correction_;}
  std::vector<D3vector>& cor_real(void){return correction_real_;}
  std::vector<D3tensor>& quad(void){return quad_;}

  // compute cos and sin matrices
  void update(void);

  ElectricEnthalpy(Sample& s);
  ~ElectricEnthalpy(void);
};
#endif
