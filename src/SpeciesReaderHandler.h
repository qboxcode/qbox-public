////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReaderHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesReaderHandler.h,v 1.1 2003-03-27 22:05:59 fgygi Exp $

#ifndef SPECIESREADERHANDLER_H
#define SPECIESREADERHANDLER_H

#include "Species.h"

#include <xercesc/sax2/DefaultHandler.hpp>
#include "StrX.h"
using namespace xercesc;
#include <sstream>
#include <vector>
#include <string>
using namespace std;

class SpeciesReaderHandler : public DefaultHandler
{
  public:
    SpeciesReaderHandler(Species& sp);
    ~SpeciesReaderHandler();

    // -----------------------------------------------------------------------
    //  Implementations of the SAX DocumentHandler interface
    // -----------------------------------------------------------------------
    void endDocument();

    void endElement(const XMLCh* const uri, const XMLCh* const localname, 
                    const XMLCh* const qname);
    void characters(const XMLCh* const chars, const unsigned int length);
    void ignorableWhitespace(const XMLCh* const chars, 
      const unsigned int length);
    void processingInstruction(const XMLCh* const target,
         const XMLCh* const data);
    void startDocument();
    void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes);

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

    int current_l, current_size;
    double current_interval;
    
  private :

    string st;
    Species& sp_;
    
};
#endif
