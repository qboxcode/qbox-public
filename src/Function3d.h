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

  D3vector a, b, c;  // domain size
  int nx, ny, nz;    // grid size
  std::string name;
  std::vector<double> val;
  void read(std::string filename);
  void write(std::string filename) const;
  void print(std::ostream& os) const;
};
std::ostream& operator << ( std::ostream& os, const Function3d& f );
#endif
