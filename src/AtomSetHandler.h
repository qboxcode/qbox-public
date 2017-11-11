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
// AtomSetHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ATOMSETHANDLER_H
#define ATOMSETHANDLER_H

#include "StructureHandler.h"
#include "D3vector.h"
#include <string>

class AtomSet;
class Species;

class AtomSetHandler : public StructureHandler
{
  private:

  AtomSet& as_;
  std::string current_atom_name, current_atom_species;
  std::string current_species_name;
  D3vector current_atom_position, current_atom_velocity;
  Species* current_species;

  public:

  // Start of the root element in the structure being handled
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of the root element in the structure being handled
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  AtomSetHandler(AtomSet& as);
  ~AtomSetHandler();
};
#endif
