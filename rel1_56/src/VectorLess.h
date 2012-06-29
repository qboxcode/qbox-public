////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// VectorLess.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VECTORLESS_H
#define VECTORLESS_H
#include <vector>

template <class T>
struct VectorLess
{
  // function object for indirect comparison of vector elements
  public:
  std::vector<T>& a_;
  VectorLess<T>(std::vector<T>& a) : a_(a) {};
  bool operator() (int i, int j) const
  {
    return a_[i] < a_[j];
  }
};
#endif
