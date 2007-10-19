////////////////////////////////////////////////////////////////////////////////
//
// SampleReader.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleReader.h,v 1.5 2007-10-19 17:37:06 fgygi Exp $

#ifndef SAMPLEREADER_H
#define SAMPLEREADER_H

enum event_type { species, atom, wavefunction, wavefunction_velocity,
                  slater_determinant, end, invalid };

class Context;
class Sample;

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
