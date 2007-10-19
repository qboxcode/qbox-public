////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.h,v 1.5 2007-10-19 16:24:05 fgygi Exp $

#ifndef SPECIESREADER_H
#define SPECIESREADER_H

#include <string>
#include "Context.h"

class SpeciesReader
{
  private:

  const Context& ctxt_;

  std::string uri_;   // uri from which Species is read

  public:

  SpeciesReader(const Context& ctxt);
  void readSpecies(Species& sp, const std::string uri);
  void bcastSpecies(Species& sp);
};

class SpeciesReaderException {};

#endif
