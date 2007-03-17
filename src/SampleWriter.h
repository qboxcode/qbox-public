////////////////////////////////////////////////////////////////////////////////
//
// SampleWriter.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleWriter.h,v 1.2 2007-03-17 01:14:00 fgygi Exp $

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
  void writeSample(const Sample& s, const std::string filename,
                   std::string description,
                   bool base64, bool atomsonly);
};

class SampleWriterException
{
  public:
  std::string msg;
  SampleWriterException(std::string s) : msg(s) {}
};

#endif
