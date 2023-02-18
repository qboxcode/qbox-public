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

typedef unsigned char byte_t;

class Base64Transcoder
{
  char etable[64];  // encode table
  byte_t dtable[256]; // decode table

  public:

  Base64Transcoder();
  void encode(size_t nbytes, const byte_t* const from, char* const to);
  size_t decode(size_t nchars, const char* const from, byte_t* const to);
  void byteswap_double(size_t n, double* const x);
  void byteswap_int(size_t n, int* const x);
  void print(size_t nchars, const char* const buf, std::ostream& o);
  void print(const std::string buf, std::ostream& o);
  void print(size_t nchars, const char* const buf, FILE* outfile);
  void print(const std::string buf, FILE* outfile);

  // number of chars needed to encode nbytes bytes
  size_t nchars(size_t nbytes) { return 4 * ( ( nbytes + 2 ) / 3 ); }
  // number of bytes needed to decode nchars chars
  size_t nbytes(size_t nchars) { return 3 * ( nchars / 4 ); }
};

#endif
