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
// MLWFTransform.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MLWFTransform.h,v 1.4 2008-08-13 06:39:43 fgygi Exp $

#ifndef MLWFTRANSFORM_H
#define MLWFTRANSFORM_H

#include <vector>
#include <complex>
class SlaterDet;
class UnitCell;
class DoubleMatrix;
#include "D3vector.h"
#include "BasisMapping.h"

class MLWFTransform
{
  private:

  const SlaterDet& sd_;
  const UnitCell& cell_;
  const Context& ctxt_;

  BasisMapping bm_;
  std::vector<DoubleMatrix*> a_;  // cosine and sine matrices
  DoubleMatrix* u_;               // orthogonal transformation
  std::vector<std::vector<double> > adiag_; // diagonal elements

  public:

  void compute_transform(void);
  void compute_sincos(const int n, const std::complex<double>* f,
    std::complex<double>* fc, std::complex<double>* fs);
  void apply_transform(SlaterDet& sd);

  double spread2(int i, int j);
  double spread2(int i);
  double spread2(void);
  double spread(int i);
  double spread(void);
  D3vector center(int i);
  D3vector dipole(void);

  MLWFTransform(const SlaterDet& sd);
  ~MLWFTransform(void);
};
#endif
