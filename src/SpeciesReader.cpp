////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2012 The Regents of the University of California
//
// This file is part of Qbox //
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SpeciesReader.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Species.h"
#include "SpeciesReader.h"
#include <cassert>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "StructuredDocumentHandler.h"
#include "SpeciesHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SpeciesReader::SpeciesReader(void) {}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::uri_to_species(const string uri, Species& sp)
{
  // parse the document defined by its URI and initialize
  // the Species object sp

  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    cout << "  Species::readSpecies: Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << endl;
    throw;
  }

  parser = XMLReaderFactory::createXMLReader();
  if (valScheme == SAX2XMLReader::Val_Auto)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, true);
  }

  if (valScheme == SAX2XMLReader::Val_Never)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
  }

  if (valScheme == SAX2XMLReader::Val_Always)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
  }

  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
  parser->setFeature(XMLUni::fgXercesSchema, doSchema);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
  parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

  int nlink = 0;
  string current_uri = uri;

  while ( sp.description() == "undefined" && nlink++ < 5 )
  {

    try
    {
      SpeciesHandler* sp_handler = new SpeciesHandler(sp);
      StructuredDocumentHandler handler(sp_handler);
      parser->setContentHandler(&handler);
      parser->setErrorHandler(&handler);
      parser->parse(uri.c_str());
      // errorCount = parser->getErrorCount();
      delete sp_handler;
    }

    catch (const XMLException& toCatch)
    {
      cout << "\nAn error occurred\n  Error: "
           << StrX(toCatch.getMessage())
           << "\n" << endl;
      XMLPlatformUtils::Terminate();
      delete parser;
      throw;
    }

    catch (const SAXParseException& e)
    {
      cout << "\na SAXParseException occurred in file "
           << StrX(e.getSystemId())
           << ", line " << e.getLineNumber()
           << ", char " << e.getColumnNumber()
           << "\n  Message: " << StrX(e.getMessage()) << endl;
      XMLPlatformUtils::Terminate();
      delete parser;
      throw;
    }

    if ( sp.description() == "undefined" )
      current_uri = sp.uri();
  }
  sp.uri_ = current_uri;

  delete parser;
  XMLPlatformUtils::Terminate();
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::uri_to_string(const string uri, const string name,
     string& xmlstr)
{
  // parse a species document defined by its URI and initialize
  // the string xmlstr with the corresponding XML <species> element
  //
  Species sp(name);
  uri_to_species(uri,sp);
  ostringstream spstream;
  const bool expanded_form = true;
  sp.print(spstream,expanded_form);
  xmlstr = spstream.str();
}

////////////////////////////////////////////////////////////////////////////////
void SpeciesReader::string_to_species(const string xmlstr, Species& sp)
{
  // parse the xmlstr string containing a <species> element and
  // initialize the Species object sp
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;

  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    cout << "  Species::readSpecies: Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << endl;
    throw;
  }

  parser = XMLReaderFactory::createXMLReader();
  if (valScheme == SAX2XMLReader::Val_Auto)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, true);
  }

  if (valScheme == SAX2XMLReader::Val_Never)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
  }

  if (valScheme == SAX2XMLReader::Val_Always)
  {
    parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
    parser->setFeature(XMLUni::fgXercesDynamic, false);
  }

  parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
  parser->setFeature(XMLUni::fgXercesSchema, doSchema);
  parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
  parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

  MemBufInputSource* memBufIS = 0;

  try
  {
    memBufIS = new MemBufInputSource
    ( (const XMLByte*) &xmlstr[0], xmlstr.size(), "buf_id", false);
    SpeciesHandler* sp_handler = new SpeciesHandler(sp);
    StructuredDocumentHandler handler(sp_handler);
    parser->setContentHandler(&handler);
    parser->setErrorHandler(&handler);
    parser->parse(*memBufIS);
    // errorCount = parser->getErrorCount();
    delete sp_handler;
  }

  catch (const XMLException& toCatch)
  {
    cout << "\nAn error occurred\n  Error: "
         << StrX(toCatch.getMessage())
         << "\n" << endl;
    XMLPlatformUtils::Terminate();
    delete parser;
    throw;
  }

  catch (const SAXParseException& e)
  {
    cout << "\na SAXParseException occurred in file "
         << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
    XMLPlatformUtils::Terminate();
    delete parser;
    throw;
  }

  delete parser;
  XMLPlatformUtils::Terminate();
}
