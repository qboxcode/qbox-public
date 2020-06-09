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
// SpeciesHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SpeciesHandler.h"
#include "Species.h"
#include "StrX.h"
using namespace xercesc;
#include <iostream>
#include <sstream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SpeciesHandler::SpeciesHandler(Species& sp) :
  sp_(sp), d_ij_alloc(false) {}

////////////////////////////////////////////////////////////////////////////////
SpeciesHandler::~SpeciesHandler() {}

////////////////////////////////////////////////////////////////////////////////
// read attributes
void SpeciesHandler::read(const Attributes& attributes)
{
  unsigned int len = attributes.getLength();
  for ( unsigned int index = 0; index < len; index++ )
  {
    string attrname = StrX(attributes.getLocalName(index)).localForm();
    if ( attrname == "l" )
    {
      current_l = atoi(StrX(attributes.getValue(index)).localForm());
    }
    else if ( attrname == "size" )
    {
      current_size = atoi(StrX(attributes.getValue(index)).localForm());
    }
    else if ( attrname == "i" )
    {
      current_i = atoi(StrX(attributes.getValue(index)).localForm()) - 1;
    }
    else if ( attrname == "j" )
    {
      current_j = atoi(StrX(attributes.getValue(index)).localForm()) - 1;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// allocate array for d_ij matrix
void SpeciesHandler::alloc_d_ij()
{
  // allocate for every l
  sp_.d_.resize(sp_.proj_.size());
  for ( int l = 0; l < sp_.d_.size(); l++ )
  {
    const int size = sp_.proj_[l].size();
    sp_.d_[l].resize(size);
    for ( int i = 0; i < size; i++ )
    {
      sp_.d_[l][i].resize(size);
    }
  }
  d_ij_alloc = true;
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
//  cout << " SpeciesHandler::startElement " << StrX(qname) << endl;

  string locname = StrX(localname).localForm();

  if ( locname == "species" )
  {
    // check for the case where the species is a link to another uri
    unsigned int len = attributes.getLength();
    for ( unsigned int index = 0; index < len; index++ )
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "name" )
      {
        current_name = StrX(attributes.getValue(index)).localForm();
        sp_.name_ = current_name;
      }
      else if ( attrname == "href" )
      {
        current_href = StrX(attributes.getValue(index)).localForm();
        sp_.uri_ = current_href;
      }
    }
  }
  else if ( locname == "norm_conserving_pseudopotential" )
  {
    sp_.type_ = Species::NCPP;
  }
  else if ( locname == "norm_conserving_semilocal_pseudopotential" )
  {
    sp_.type_ = Species::SLPP;
  }
  else if ( locname == "core_density" )
  {
    read(attributes);
  }
  else if ( locname == "local_potential" )
  {
    read(attributes);
  }
  else if ( locname == "projector" )
  {
    read(attributes);
  }
  else if ( locname == "d_ij" )
  {
    read(attributes);
  }
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  string locname = StrX(localname).localForm();
  istringstream stst(content);

  if ( locname == "description" )
  {
    // reject ambiguous case where both the href and the definition are given
    if ( current_href != "" )
    {
      cout << " SpeciesHandler: ambiguous definition: uri=" << StrX(uri)
        << endl << " using local definition (href: " << current_href
        << " ignored)" << endl;
    }
    sp_.description_ = content;
  }
  else if ( locname == "atomic_number" )
  {
    stst >> sp_.atomic_number_;
  }
  else if ( locname == "mass" )
  {
    stst >> sp_.mass_;
  }
  else if ( locname == "symbol" )
  {
    stst >> skipws >> sp_.symbol_;
  }
  else if ( locname == "valence_charge" )
  {
    stst >> sp_.zval_;
  }
  else if ( locname == "lmax" )
  {
    stst >> sp_.lmax_;
  }
  else if ( locname == "llocal" )
  {
    stst >> sp_.llocal_;
  }
  else if ( locname == "nquad" )
  {
    stst >> sp_.nquad_;
  }
  else if ( locname == "rquad" )
  {
    stst >> sp_.rquad_;
  }
  else if ( locname == "mesh_spacing" )
  {
    stst >> sp_.deltar_;
  }
  else if ( locname == "radial_potential" )
  {
    assert(sp_.type_ == Species::NCPP);
    if ( current_l + 1 > sp_.vps_.size() )
    {
      sp_.vps_.resize(current_l + 1);
      sp_.phi_.resize(current_l + 1);
    }
    sp_.vps_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.vps_[current_l][i];
  }
  else if ( locname == "local_potential" )
  {
    assert(sp_.type_ == Species::SLPP);
    sp_.vlocal_.resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.vlocal_[i];
  }
  else if ( locname == "radial_function" )
  {
    assert(sp_.type_ == Species::NCPP);
    sp_.phi_[current_l].resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.phi_[current_l][i];
  }
  else if ( locname == "core_density" )
  {
    // read charge for nonlinear core correction
    sp_.nlcc_.resize(current_size);
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.nlcc_[i];
  }
  else if ( locname == "projector" and sp_.type_ == Species::SLPP )
  {
    // read one of the projector with this angular momentum
    // resize vector if necessary
    if ( current_l >= sp_.proj_.size() ) sp_.proj_.resize(current_l + 1);
    if ( current_i >= sp_.proj_[current_l].size() )
      sp_.proj_[current_l].resize(current_i + 1);
    sp_.proj_[current_l][current_i].resize(current_size);
    // read the projector
    for ( int i = 0; i < current_size; i++ )
      stst >> sp_.proj_[current_l][current_i][i];
  }
  else if ( locname == "d_ij" )
  {
    assert(sp_.type_ == Species::SLPP);
    if ( not d_ij_alloc ) alloc_d_ij();
    assert(current_l < sp_.d_.size());
    assert(current_i < sp_.d_[current_l].size());
    assert(current_j < sp_.d_[current_l][current_i].size());
    stst >> sp_.d_[current_l][current_i][current_j];
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* SpeciesHandler::startSubHandler(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " SpeciesHandler::startSubHandler " << StrX(qname) << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesHandler::endSubHandler(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const StructureHandler* const subHandler)
{
  // cout << " SpeciesHandler::endSubHandler " << StrX(qname) << endl;
  // if any StructureHandler was created by startSubHandler, delete it
  // delete subHandler;
}
