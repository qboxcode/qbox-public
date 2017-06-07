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
// NameOf.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NAMEOF_H
#define NAMEOF_H

#include <string>

// predicate class for searching T* containers by name
// T must be a pointer type to something that has a name() member
template <class T>
class NameOf
{
  public:
  std::string name;

  NameOf<T>(std::string s) : name(s) {};
  bool operator() (T t) const
  {
    return t->name() == name;
  }
};
#endif
