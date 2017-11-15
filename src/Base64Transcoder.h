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
// Base64Transcoder.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BASE64TRANSCODER_H
#define BASE64TRANSCODER_H

#include <iostream>
#include <cstdio>
#include <string>

typedef unsigned char byte;

class Base64Transcoder
{
  char etable[64];  // encode table
  byte dtable[256]; // decode table

  public:

  Base64Transcoder();
  int encode(int nbytes, const byte* const from, char* const to);
  int decode(int nchars, const char* const from, byte* const to);
  void byteswap_double(size_t n, double* const x);
  void byteswap_int(size_t n, int* const x);
  int print(int nchars, const char* const buf, std::ostream& o);
  int print(const std::string buf, std::ostream& o);
  int print(int nchars, const char* const buf, FILE* outfile);
  int print(const std::string buf, FILE* outfile);

  // number of chars needed to encode nbytes bytes
  int nchars(int nbytes) { return 4 * ( ( nbytes + 2 ) / 3 ); }
  // number of bytes needed to decode nchars chars
  int nbytes(int nchars) { return 3 * ( nchars / 4 ); }
};

#endif
