////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2012 The Regents of the University of California
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

#ifndef SPECIESREADER_H
#define SPECIESREADER_H

#include <string>

class SpeciesReader
{
  private:

  public:

  SpeciesReader(void);

  // initialize a Species object using a species URI
  void uri_to_species(const std::string uri, Species& sp);

  // initialize a string containing an XML <species> element using a species URI
  void uri_to_string(const std::string uri, const std::string name,
       std::string& xmlstr);

  // initialize a Species object using a string containing a <species> element
  void string_to_species(const std::string xmlstr, Species& sp);
};
#endif
