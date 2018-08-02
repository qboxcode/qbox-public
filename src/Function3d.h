////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 The Regents of the University of California
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
// Function3d.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FUNCTION3D_H
#define FUNCTION3D_H

#include "D3vector.h"
#include <vector>
#include <string>
#include <map>

class Function3d
{
  public:

  D3vector a_, b_, c_;
  int gnx_, gny_, gnz_, nx_, ny_, nz_, x0_, y0_, z0_;
  std::string str_;
  std::string name_;
  std::vector<double> val_;

  const D3vector a(void) const { return a_; }
  const D3vector b(void) const { return b_; }
  const D3vector c(void) const { return c_; }
  int nx(void) const { return nx_; }
  int ny(void) const { return ny_; }
  int nz(void) const { return nz_; }
  int size(void) const { return val_.size(); }

  std::vector<double>& val(void) { return val_; }

  void read(std::string filename);
  void write(std::string filename, std::string encoding) const;
};
#endif
