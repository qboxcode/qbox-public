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
// SampleReader.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "SampleReader.h"

#include "XMLGFPreprocessor.h"

#include "Timer.h"

#include <cassert>
#include <string>
#include <cstring> // memcpy
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "SampleHandler.h"
#include "StructuredDocumentHandler.h"
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void SampleReader::readSample (Sample& s, const string uri, bool serial)
{
  Timer tm;
  tm.start();
  // create a global Context with MPIdata::ngb() rows
  Context gctxt(MPI_COMM_WORLD,MPIdata::ngb(),MPIdata::size()/MPIdata::ngb());
  const bool onpe0 = gctxt.onpe0();

  SAX2XMLReader::ValSchemes valScheme;
  valScheme = SAX2XMLReader::Val_Auto;
  //alternate choices of validation schemes
  //valScheme = SAX2XMLReader::Val_Always;
  //valScheme = SAX2XMLReader::Val_Never;

  bool doNamespaces = true;
  // Validate agains schema on pe0 only to avoid multiple warning messages
  // This does not weaken validation since all tasks parse the same document
  bool doSchema = onpe0;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  MemBufInputSource* memBufIS = 0;

  // XMLPlatformUtils initialization on all tasks
  try
  {
    XMLPlatformUtils::Initialize();
  }

  catch (const XMLException& toCatch)
  {
    cout << "  SampleReader::readSample: Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << endl;
    return;
  }

  string xmlcontent;
  DoubleMatrix gfdata(gctxt);

  XMLGFPreprocessor xmlgfp;
  if ( xmlgfp.process(uri.c_str(),gfdata,xmlcontent,serial) )
  {
    cout << "  SampleReader::readSample: Error in XMLGFPreprocessor" << endl;
    return;
  }

  // Each task holds a copy of xmlcontent
  // The distributed matrix gfdata contains the grid function data

  if ( onpe0 )
  {
    cout << " xmlcontent.size(): " << xmlcontent.size() << endl;
#if DEBUG
    cout << MPIdata::rank() << ": xmlcontent: " << xmlcontent << endl;
#endif
  }
  memBufIS = new MemBufInputSource
    ( (const XMLByte*) &xmlcontent[0], xmlcontent.size(), "buf_id", false);

  // parse xmlcontent on each task
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

  SampleHandler* s_handler = new SampleHandler(s,gfdata);

  try
  {
    StructuredDocumentHandler handler(s_handler);
    parser->setContentHandler(&handler);
    parser->setErrorHandler(&handler);

    gctxt.barrier();
    if ( onpe0 ) cout << " Starting XML parsing" << endl;
      parser->parse(*memBufIS);
    if ( onpe0 ) cout << " XML parsing done" << endl;
    // errorCount = parser->getErrorCount();
  }

  catch (const XMLException& toCatch)
  {
    cout << "\nAn error occurred\n  Error: "
         << StrX(toCatch.getMessage())
         << "\n" << endl;
  }

  catch (const SAXParseException& e)
  {
    cout << "\na SAXParseException occurred in file "
         << StrX(e.getSystemId())
         << ", line " << e.getLineNumber()
         << ", char " << e.getColumnNumber()
         << "\n  Message: " << StrX(e.getMessage()) << endl;
  }

  delete s_handler;
  delete parser;
  delete memBufIS;
  XMLPlatformUtils::Terminate();

  tm.stop();
  if ( onpe0 )
    cout << " SampleReader: read time: " << tm.real() << " s" << endl;

  // Enforce consistency between Wavefunction and AtomSet

  // If a Wavefunction was read, or if the Wavefunction was
  // previously defined, the Wavefunction cell volume is non-zero.
  // In this case, enforce consistency between Wavefunction and
  // AtomSet cell and nel
  if ( s.wf.cell().volume() != 0.0 )
  {
    // copy the Wavefunction cell to the AtomSet cell
    s.atoms.set_cell(s.wf.cell());

    // if the number of electrons differ between the new AtomSet and
    // the Wavefunction, throw an exception
    if ( s.atoms.nel() != s.wf.nel() )
      throw runtime_error("number of electrons incompatible");
  }

  // If the Wavefunction cell volume is zero, the Wavefunction was
  // not read and was not previously set,
  // In this case, initialize the Wavefunction using the values of
  // cell and nel from the AtomSet
  if ( s.wf.cell().volume() == 0.0 )
  {
    // Initialize the Wavefunction
    s.wf.reset();
    // set wf cell
    s.wf.resize(s.atoms.cell(),s.atoms.cell(),s.wf.ecut());
    // set number of states from charge in atomset
    s.wf.set_nel(s.atoms.nel());
    s.wf.update_occ(0.0);

    // if a wavefunction_velocity was read, delete it
    if ( s.wfv != 0 )
    {
      delete s.wfv;
      s.wfv = 0;
    }
  }
}
