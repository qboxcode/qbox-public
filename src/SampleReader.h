////////////////////////////////////////////////////////////////////////////////
//
// SampleReader.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleReader.h,v 1.4 2007-10-19 16:24:05 fgygi Exp $

#ifndef SAMPLEREADER_H
#define SAMPLEREADER_H

class Context;

class SampleReader
{
  private:

  const Context& ctxt_;

  public:

  SampleReader(const Context& ctxt);
  void readSample(Sample& s, const std::string uri, bool serial);
};

class SampleReaderException
{
  public:
  std::string msg;
  SampleReaderException(std::string s) : msg(s) {}
};

#endif
