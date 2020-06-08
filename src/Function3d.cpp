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
#include "Function3dHandler.h"
#include "Base64Transcoder.h"
#include "qbox_xmlns.h"
#include "Timer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstdio>
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
  cout << "Function3d::read: file size: " << sz << endl;

  // read file into xmlcontent string
  string xmlcontent;
  xmlcontent.resize(sz);

  tm.start();
  size_t items_read = fread(&xmlcontent[0],sizeof(char),sz,infile);
  assert(items_read==sz);
  fclose(infile);
  tm.stop();

  cout << "Function3d::read: fread time: " << tm.real() << endl;
  cout << "Function3d::read: read rate: "
       << sz/(tm.real()*1024*1024) << " MB/s" << endl;

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
  cout << "Function3d::read: grid size: "
       << nx << " " << ny << " " << nz << endl;
  cout << "Function3d::read: val size: " << val.size() << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Function3d::write(string filename) const
{
  ofstream os(filename.c_str());
  print(os);
}

////////////////////////////////////////////////////////////////////////////////
void Function3d::print(ostream& os) const
{
  Base64Transcoder xcdr;
  int nbytes = val.size() * sizeof(double);
  int nchars = xcdr.nchars(nbytes);
  char *wbuf = new char[nchars];
#if PLT_BIG_ENDIAN
  // make copy of val for byte swapping without affecting original array
  vector<double> tmpval(val);
  xcdr.byteswap_double(tmpval.size(),&tmpval[0]);
  // transform tmpval to base64 encoding
  xcdr.encode(nbytes, (byte *) &tmpval[0], wbuf);
#else
  // transform val to base64 encoding
  xcdr.encode(nbytes, (byte *) &val[0], wbuf);
#endif

  os <<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
     <<"<fpmd:function3d xmlns:fpmd=\""
     << qbox_xmlns()
     << "\"\n"
     <<" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
     <<" xsi:schemaLocation=\""
     << qbox_xmlns() << " function3d.xsd\"\n"
     << " name=\"" << name << "\">"
     << endl;
  os.setf(ios::fixed,ios::floatfield);
  os.precision(8);
  os << "<domain a=\"" << a << "\"" << endl;
  os << "        b=\"" << b << "\"" << endl;
  os << "        c=\"" << c << "\"/>" << endl;
  os << "<grid nx=\"" << nx << "\" ny=\"" << ny
     << "\" nz=\"" << nz << "\"/>" << endl;
  os << "<grid_function type=\"double\" nx=\"" << nx
     << "\" ny=\"" << ny << "\" nz=\"" << nz
     << "\" encoding=\"base64\">" << endl;
  xcdr.print(nchars,wbuf,os);
  os << "</grid_function>" << endl;
  os << "</fpmd:function3d>" << endl;
  delete [] wbuf;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const Function3d& f)
{
  f.print(os);
  return os;
}
