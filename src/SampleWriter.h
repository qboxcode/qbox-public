////////////////////////////////////////////////////////////////////////////////
//
// SampleWriter.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleWriter.h,v 1.1 2007-01-27 23:43:55 fgygi Exp $

#ifndef SAMPLEWRITER_H
#define SAMPLEWRITER_H

#include "Context.h"
class Sample;

class SampleWriter
{
  private:

  const Context& ctxt_;
  
  public:

  SampleWriter(const Context& ctxt);
  void writeSample(const Sample& s, const string filename,
                   string description,
                   bool base64, bool atomsonly);
};

class SampleWriterException
{
  public:
  string msg;
  SampleWriterException(string s) : msg(s) {}
};

#endif
