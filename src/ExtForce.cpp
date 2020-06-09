////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
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
// ExtForce.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "ExtForce.h"
#include <iostream>
using namespace std;

ostream& operator << ( ostream &os, ExtForce &x )
{
  return x.print(os);
}

