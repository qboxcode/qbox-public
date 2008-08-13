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
// SpeciesReader.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReader.h,v 1.6 2008-08-13 06:39:43 fgygi Exp $

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
