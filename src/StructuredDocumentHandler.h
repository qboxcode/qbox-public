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
// StructuredDocumentHandler.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STRUCTUREDDOCUMENTHANDLER_H
#define STRUCTUREDDOCUMENTHANDLER_H

#include <xercesc/sax2/DefaultHandler.hpp>
#include "StrX.h"

#include "StructureHandler.h"

#include <stack>
#include <string>

class StructuredDocumentHandler : public DefaultHandler
{
  struct HandlerContext
  {
    int depth;
    StructureHandler* handler;
    HandlerContext(StructureHandler* handler_, int depth_) :
      handler(handler_), depth(depth_) {}
  };

  protected:

  std::stack<HandlerContext> contextStack;
  int nestingDepth;
  int contextDepth;
  StructureHandler* activeHandler;
  std::string buffer;

  public:

  StructuredDocumentHandler(StructureHandler* handler) :
    activeHandler(handler), contextDepth(0), nestingDepth(0) {}

  ~StructuredDocumentHandler() {}

  // -----------------------------------------------------------------------
  //  Implementations of the SAX DocumentHandler interface
  // -----------------------------------------------------------------------
  void startDocument();
  void endDocument();

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
  void ignorableWhitespace(const XMLCh* const chars,
    const unsigned int length);
  void processingInstruction(const XMLCh* const target,
       const XMLCh* const data);

  // -----------------------------------------------------------------------
  //  Implementations of the SAX ErrorHandler interface
  // -----------------------------------------------------------------------
  void warning(const SAXParseException& exception);
  void error(const SAXParseException& exception);
  void fatalError(const SAXParseException& exception);

  // -----------------------------------------------------------------------
  //  Implementation of the SAX DTDHandler interface
  // -----------------------------------------------------------------------
  void notationDecl(const XMLCh* const name, const XMLCh* const publicId,
       const XMLCh* const systemId);

  void unparsedEntityDecl(const XMLCh* const name,
    const XMLCh* const publicId, const XMLCh* const systemId,
    const XMLCh* const notationName);

};
#endif
