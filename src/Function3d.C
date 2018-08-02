//////////////////////////////////////////////////////////////////////////////// //
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
// Function3d.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Function3d.h"
#include "Timer.h"
#include "Function3dHandler.h"
#include "Base64Transcoder.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>
using namespace std;
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;

////////////////////////////////////////////////////////////////////////////////
void Function3d::read(string filename)
{
  Timer tm;
  struct stat statbuf;
  bool filestat = !stat(filename.c_str(),&statbuf);
  if ( !filestat )
  {
    cout << "Function3d::read: could not stat file " << filename << endl;
    return;
  }

  FILE* infile;
  infile = fopen(filename.c_str(),"r");
  if ( !infile )
  {
    cout << "Function3d::read: could not open file " << filename
         << " for reading" << endl;
    return;
  }

  // determine size
  off_t sz = statbuf.st_size;

  // use contiguous read buffer, to be copied later to a string
  char *rdbuf = new char[sz];
  cout << "Function3d::read: file size: " << sz << endl;

  tm.start();
  size_t items_read = fread(rdbuf,sizeof(char),sz,infile);
  assert(items_read==sz);
  fclose(infile);
  tm.stop();

  cout << "Function3d::read: fread time: " << tm.real() << endl;
  cout << "Function3d::read: read rate: "
       << sz/(tm.real()*1024*1024) << " MB/s" << endl;

  string xmlcontent;
  xmlcontent.assign(rdbuf,sz);
  delete [] rdbuf;

  tm.reset();
  tm.start();

  SAX2XMLReader::ValSchemes valScheme;
  valScheme = SAX2XMLReader::Val_Never;
  //alternate choices of validation schemes
  //valScheme = SAX2XMLReader::Val_Always;
  //valScheme = SAX2XMLReader::Val_Never;

  bool doNamespaces = true;
  bool doSchema = false;
  bool schemaFullChecking = false;
  bool namespacePrefixes = false;

  try
  {
    XMLPlatformUtils::Initialize();
  }

  catch (const XMLException& toCatch)
  {
    cout << "Function3d::read: Error during XML initialization :\n"
         << StrX(toCatch.getMessage()) << endl;
    return;
  }

  cout << "Function3d::read: xmlcontent.size(): " << xmlcontent.size() << endl;

  MemBufInputSource* memBufIS = new MemBufInputSource
    ( (const XMLByte*) &xmlcontent[0], xmlcontent.size(), "buf_id", false);

  // parse xmlcontent
  SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
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

  Function3dHandler* f_handler = new Function3dHandler(*this);

  try
  {
    parser->setContentHandler(f_handler);
    parser->setErrorHandler(f_handler);
    cout << "Function3d::read: Starting XML parsing" << endl;
    parser->parse(*memBufIS);
    cout << "Function3d::read: XML parsing done" << endl;
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

  delete f_handler;
  delete parser;
  delete memBufIS;
  XMLPlatformUtils::Terminate();

  tm.stop();
  cout << "Function3d::read: parse time: " << tm.real() << endl;

  tm.reset();
  tm.start();
  val_.resize(gnx_*gny_*gnz_);
  cout << "Function3d::read: grid size: "
       << gnx_ << " " << gny_ << " " << gnz_ << endl;
  cout << "Function3d::read: str_ size: " << str_.size() << endl;

  Base64Transcoder xcdr;
  size_t nbytes = xcdr.decode(str_.size(),str_.data(),(byte*)&val_[0]);
  tm.stop();
  assert(nbytes==val_.size()*sizeof(double));
  cout << "Function3d::read: base64 transcode time: " << tm.real() << endl;
  str_.clear();
}
