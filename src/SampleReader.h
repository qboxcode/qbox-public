////////////////////////////////////////////////////////////////////////////////
//
// SampleReader.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleReader.h,v 1.2 2003-08-22 18:01:13 fgygi Exp $

#ifndef SAMPLEREADER_H
#define SAMPLEREADER_H

#include "Context.h"

class SampleReader
{
  private:

  const Context& ctxt_;
  
  public:

  SampleReader(const Context& ctxt);
  void readSample(Sample& s, const string uri, bool serial);
};

class SampleReaderException
{
  public:
  string msg;
  SampleReaderException(string s) : msg(s) {}
};

#endif
