////////////////////////////////////////////////////////////////////////////////
//
// Base64Transcoder.h
//
////////////////////////////////////////////////////////////////////////////////
// $ Id: $

#ifndef BASE64TRANSCODER_H
#define BASE64TRANSCODER_H

#include <iostream>
#include <string>
using namespace std;
typedef unsigned char byte;

class Base64Transcoder
{
  char etable[256]; // encode table
  byte dtable[256]; // decode table

  public:
  
  Base64Transcoder();
  int encode(int nbytes, const byte* const from, char* const to);
  int decode(int nchars, const char* const from, byte* const to);
  void byteswap_double(size_t n, double* const x);
  void byteswap_int(size_t n, int* const x);
  int print(int nchars, const char* const buf, ostream& o);
  int print(const string buf, ostream& o);

  // number of chars needed to encode nbytes bytes
  int nchars(int nbytes) { return 4 * ( ( nbytes + 2 ) / 3 ); }
  // number of bytes needed to decode nchars chars
  int nbytes(int nchars) { return 3 * ( nchars / 4 ); }
};

#endif
