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
// isodate.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "isodate.h"
#include <ctime>
std::string isodate(void)
{
  const time_t t = time(NULL);
  struct tm* tms = gmtime(&t);
  char s[32];
  const char* fmt = "%Y-%m-%dT%TZ";
  strftime(s,32,fmt,tms);
  return std::string(s);
}
