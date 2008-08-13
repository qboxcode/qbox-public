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
// SpeciesHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesHandler.h,v 1.4 2008-08-13 06:39:43 fgygi Exp $

#ifndef SPECIESHANDLER_H
#define SPECIESHANDLER_H

#include "StructureHandler.h"

class Species;

class SpeciesHandler : public StructureHandler
{
  private:

  Species& sp_;
  int current_l, current_size;
  std::string current_name, current_href;
  double current_interval;

  public:

  // Start of an element handled by SpeciesHandler
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

  // End of an element handled by SpeciesHandler
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname,
      const XMLCh* const qname, std::string& content);

  // start a subHandler if possible
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes);

  // end subHandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler);

  SpeciesHandler(Species& sp);
  ~SpeciesHandler();
};
#endif
