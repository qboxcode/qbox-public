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
  Wavefunction* dwf_;
  SlaterDet& sd_;
  const Context& ctxt_;
  const Context& vctxt_;

  bool onpe0_;
  bool finite_field_;
  Basis* vbasis_;

  enum { off, berry, mlwf, mlwf_ref, mlwf_ref_q } pol_type_;
  bool compute_quadrupole_;

  // electric field
  D3vector e_field_;

  Wavefunction* rwf_[3];

  // MLWFtransform is used to compute S matrix
  MLWFTransform* mlwft_;

  // s matrices
  ComplexMatrix* smat_[3];

  // total, ionic and electronic part of macroscopic polarization
  D3vector polarization_total_, polarization_ion_, polarization_elec_;

  // polarization energy
  double energy_;

  std::vector <D3vector> mlwfc_;
  std::vector <double> mlwfs_;
  std::vector <D3vector> correction_;
  std::vector <D3tensor> quad_;

  void compute_correction(void);

  public:

  mutable TimerMap tmap;

  D3vector e_field(void) const { return e_field_; }
  D3vector polarization(void) const { return polarization_total_; }
  D3vector polarization_ion(void) const { return polarization_ion_; }
  D3vector polarization_elec(void) const { return polarization_elec_; }

  double energy(Wavefunction& dwf, bool compute_hpsi);

  void update(void);
  void print(std::ostream& os) const;

  ElectricEnthalpy(Sample& s);
  ~ElectricEnthalpy(void);
};
std::ostream& operator << ( std::ostream& os, const ElectricEnthalpy& e );
#endif
