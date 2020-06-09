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
// StructuredDocumentHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#if USE_XERCES

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include "StructuredDocumentHandler.h"

#if TIMING
#include "Timer.h"
#endif

#include "StrX.h"
#include <iostream>
#include <cassert>
using namespace xercesc;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " StructuredDocumentHandler::startElement " << StrX(qname) << endl;
  // cout << " nestingDepth before: " << nestingDepth << endl;
  buffer = "";

  // attempt to start a subhandler for this element
  StructureHandler* next = activeHandler->startSubHandler(uri,localname,qname,
    attributes);

  // yes, this element can be processed by a subhandler
  if ( next != 0 )
  {
    contextStack.push(HandlerContext(activeHandler,contextDepth));

    activeHandler = next;
    contextDepth = nestingDepth;
  }
  activeHandler->startElement(uri,localname,qname,attributes);
  nestingDepth++;
  // cout << " nestingDepth after:  " << nestingDepth << endl;
}

////////////////////////////////////////////////////////////////////////////////
#ifndef XERCES_VERSION_MAJOR
#error "XERCES_VERSION_MAJOR not defined"
#endif
#if XERCES_VERSION_MAJOR > 2
void StructuredDocumentHandler::characters(const XMLCh* const chars,
  const XMLSize_t length)
#else
void StructuredDocumentHandler::characters(const XMLCh* const chars,
  const unsigned int length)
#endif
{
#if TIMING
  Timer tm;
  tm.start();
#endif

  char *str = XMLString::transcode(chars);
  buffer += str;
  XMLString::release(&str);

#if TIMING
  tm.stop();
  cout << " StructuredDocumentHandler::characters: time: " << tm.real() << endl;
#endif
  //cout << "length=" << length << " buffer.size=" << buffer.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname)
{
  // cout << " StructuredDocumentHandler::endElement " << StrX(qname) << endl;
  // cout << " nestingDepth before: " << nestingDepth << endl;
  nestingDepth--;
  activeHandler->endElement(uri,localname,qname,buffer);

  // Check if this element was processed by a subhandler
  // If yes, pop the previous context from the Context stack
  if ( nestingDepth != 0 && nestingDepth == contextDepth )
  {
    // cout << " popping context stack: nestingDepth=" << nestingDepth << endl;
    HandlerContext context = contextStack.top();
    StructureHandler* last = activeHandler;
    activeHandler = context.handler;
    contextDepth = context.depth;
    contextStack.pop();

    // notify activeHandler that current structure is complete
    // This allows the activeHandler to delete the last handler it has created
    activeHandler->endSubHandler(uri,localname,qname,last);
  }
  // cout << " nestingDepth after:  " << nestingDepth << endl;
  buffer = "";
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::startDocument() {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::endDocument() {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::ignorableWhitespace( const   XMLCh* const chars,
  const  unsigned int length) {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::processingInstruction(const  XMLCh* const target,
  const XMLCh* const data) {}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::error(const SAXParseException& e)
{
    cout << "\nError at file " << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;

    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::fatalError(const SAXParseException& e)
{
    cout << "\nFatal Error at file " << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
    throw(e);
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::warning(const SAXParseException& e)
{
    cout << "\nWarning at file " << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::unparsedEntityDecl(const     XMLCh* const name,
  const   XMLCh* const publicId, const   XMLCh* const systemId,
  const   XMLCh* const notationName)
{
    // Not used at this time
}

////////////////////////////////////////////////////////////////////////////////
void StructuredDocumentHandler::notationDecl(const   XMLCh* const name,
  const XMLCh* const publicId, const XMLCh* const systemId)
{
    // Not used at this time
}

#endif
