////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.h,v 1.3 2003-05-23 21:51:04 fgygi Exp $

#ifndef SPECIESREADER_H
#define SPECIESREADER_H

#include <string>
using namespace std;
#include "Context.h"

class SpeciesReader
{
  private:
  
  const Context& ctxt_;
  
  string uri_;   // uri from which Species is read
  
  public:

  SpeciesReader(const Context& ctxt);
  void readSpecies(Species& sp, const string uri);
  void bcastSpecies(Species& sp);
};

class SpeciesReaderException {};

#endif
