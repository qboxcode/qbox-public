////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 The Regents of the University of California
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
// Function3dHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FUNCTION3DHANDLER_H
#define FUNCTION3DHANDLER_H

#include "Function3d.h"
#include <xercesc/sax2/DefaultHandler.hpp>
#include "StrX.h"
#include <string>

class Function3dHandler: public DefaultHandler
{
  private:

  std::string buf_;
  Function3d& f_;
  int fnx_, fny_, fnz_; // size of fragment read
  int x0_, y0_, z0_;    // offset of fragment read

  public:

  Function3dHandler(Function3d& f);
  virtual ~Function3dHandler(void);

  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------

  void startElement(const XMLCh* const uri,const XMLCh* const localname,
    const XMLCh* const qname, const Attributes& attributes);
#ifndef XERCES_VERSION_MAJOR
#error "XERCES_VERSION_MAJOR not defined"
#endif
#if XERCES_VERSION_MAJOR > 2
  void characters(const XMLCh* const chars, const XMLSize_t length);
#else
  void characters(const XMLCh* const chars, const unsigned int length);
#endif
  void endElement(const XMLCh* const uri, const XMLCh* const localname,
                  const XMLCh* const qname);
};
#endif
