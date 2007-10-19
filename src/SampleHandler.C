////////////////////////////////////////////////////////////////////////////////
//
// SampleHandler.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleHandler.C,v 1.7 2007-10-19 16:24:04 fgygi Exp $

#if USE_XERCES

#include "SampleHandler.h"
#include "Sample.h"
#include "AtomSetHandler.h"
#include "WavefunctionHandler.h"
#include "StrX.h"
using namespace xercesc;
#include <iostream>
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleHandler::SampleHandler(Sample& s, DoubleMatrix& gfdata,
  Wavefunction& wfvtmp) :
  s_(s), gfdata_(gfdata), read_wf(false), read_wfv(false), wfvtmp_(wfvtmp) {}

////////////////////////////////////////////////////////////////////////////////
SampleHandler::~SampleHandler() {}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  // cout << " SampleHandler::startElement " << StrX(qname) << endl;
}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  //istringstream stst(st);
  string locname(XMLString::transcode(localname));
  // cout << " SampleHandler::endElement " << locname << endl;
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* SampleHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " SampleHandler::startSubHandler " << StrX(qname) << endl;

  string qnm = XMLString::transcode(qname);
  if ( qnm == "atomset" )
  {
    return new AtomSetHandler(s_.atoms);
  }
  else if ( qnm == "wavefunction" )
  {
    read_wf = true;
    return new WavefunctionHandler(s_.wf,gfdata_);
  }
  else if ( qnm == "wavefunction_velocity" )
  {
    read_wfv = true;
    return new WavefunctionHandler(wfvtmp_,gfdata_);
  }
  else
  {
    return 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void SampleHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const subHandler)
{
  // cout << " SampleHandler::endSubHandler " << StrX(qname) << endl;
  delete subHandler;
}
#endif
