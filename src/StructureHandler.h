////////////////////////////////////////////////////////////////////////////////
//
// StructureHandler.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: StructureHandler.h,v 1.1 2003-05-16 16:14:00 fgygi Exp $

#ifndef STRUCTUREHANDLER_H
#define STRUCTUREHANDLER_H

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
using namespace xercesc;

#include <string>
using namespace std;

class StructureHandler
{
  private:

  public:
  
  // Start of an element handled by the StructureHandler
  virtual void startElement(const XMLCh* const uri,const XMLCh* const localname,
      const XMLCh* const qname, const Attributes& attributes) = 0;

  // End of an element handled by the StructureHandler
  virtual void endElement(const XMLCh* const uri, const XMLCh* const localname, 
      const XMLCh* const qname, string& content) = 0;
  
  // start a subhandler
  virtual StructureHandler* startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname, 
    const Attributes& attributes) = 0;
    
  // end a subhandler
  virtual void endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname, 
    const StructureHandler* const subHandler) = 0;
      
  StructureHandler() {}
  virtual ~StructureHandler() {}
};
#endif
