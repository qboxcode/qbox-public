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
// AtomSetHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "AtomSetHandler.h"
#include "AtomSet.h"
#include "Species.h"
#include "SpeciesHandler.h"
#include "SpeciesReader.h"
#include "StrX.h"
#include "SampleReader.h"
using namespace xercesc;
#include <iostream>
#include <cassert>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
AtomSetHandler::AtomSetHandler(AtomSet& as) :
  as_(as) {}

////////////////////////////////////////////////////////////////////////////////
AtomSetHandler::~AtomSetHandler() {}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " AtomSetHandler::startElement " << StrX(qname) << endl;
  string locname = StrX(localname).localForm();

  // consider only elements that are dealt with directly by AtomSetHandler
  // i.e. "atom". The "species" element is delegated to a SpeciesHandler
  if ( locname == "unit_cell")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "a")
      {
        stst >> a;
      }
      else if ( attrname == "b" )
      {
        stst >> b;
      }
      else if ( attrname == "c" )
      {
        stst >> c;
      }
    }

    as_.set_cell(a,b,c);
  }
  else if ( locname == "atom")
  {
    // set default velocity to zero
    current_atom_velocity = D3vector(0.0,0.0,0.0);

    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "name")
      {
        current_atom_name = StrX(attributes.getValue(index)).localForm();
      }
      else if ( attrname == "species" )
      {
        current_atom_species = StrX(attributes.getValue(index)).localForm();
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  string locname = StrX(localname).localForm();
  // cout << " AtomSetHandler::endElement " << locname << endl;
  istringstream stst(content);
  if ( locname == "unit_cell")
  {
  }
  else if ( locname == "atom")
  {
    // create an instance of Atom using the current values read so far
    // add the Atom to the current AtomSet
    // cout << " AtomSetHandler::endElement: creating Atom(position="
    //      << current_atom_position
    //      << ", velocity=" << current_atom_velocity
    //      << ", name=" << current_atom_name
    //      << ", species=" << current_atom_species
    //      << ")" << endl;

    try
    {
      Atom* a = new Atom(current_atom_name, current_atom_species,
                         current_atom_position, current_atom_velocity);
      as_.addAtom(a);
    }
    catch (...)
    {
      cout << " AtomSetHandler::endElement: and exception occurred"
           << " in AtomSet::addAtom" << endl;
      throw;
    }

  }
  else if ( locname == "position" )
  {
    stst >> current_atom_position;
  }
  else if ( locname == "velocity" )
  {
    stst >> current_atom_velocity;
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* AtomSetHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " AtomSetHandler::startSubHandler " << StrX(qname) << endl;

  string locname = StrX(localname).localForm();
  if ( locname == "species")
  {
    // check for species attributes
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "name")
      {
        current_species_name = StrX(attributes.getValue(index)).localForm();
      }
    }

    // delegate to SpeciesHandler
    current_species = new Species(current_species_name);
    return new SpeciesHandler(*current_species);
  }
  else
  {
    return 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void AtomSetHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const last)
{
  string locname = StrX(localname).localForm();
  // cout << " AtomSetHandler::endSubHandler " << locname << endl;
  if ( locname == "species" )
  {
    SpeciesReader sp_reader;

    // check if only the uri was provided
    if ( current_species->uri() != "" )
    {
      // href was found in species definition
      // attempt to read the species from that uri

      sp_reader.uri_to_species(current_species->uri(),*current_species);
    }

    // cout << "AtomSetHandler::endSubHandler: adding Species:"
    //      << current_species_name << endl;
    try
    {
      as_.addSpecies(current_species,current_species_name);
    }
    catch (...)
    {
      cout << " AtomSetHandler::endSubHandler: and exception occurred"
           << " in AtomSet::addSpecies" << endl;
      throw;
    }
  }
  delete last;
}
