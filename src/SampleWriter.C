////////////////////////////////////////////////////////////////////////////////
//
// SampleWriter.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleWriter.C,v 1.1 2007-01-27 23:43:55 fgygi Exp $


#include "SampleWriter.h"
#include "Sample.h"
#include "fstream"
#include "qbox_xmlns.h"

#ifdef USE_CSTDIO_LFS
#include <cstdio>
#include <cstdlib>
#include <sstream>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleWriter::SampleWriter(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SampleWriter::writeSample(const Sample& s, const string filename,
                              string description,
                              bool base64, bool atomsonly)
{
  // set default encoding
  string encoding =  base64 ? "base64" : "text";
  const char* filename_cstr = filename.c_str();
  
#if USE_CSTDIO_LFS
  // This section for compilers with broken large file support
  // As of 2003-09-30, this includes gcc-3.2, icc-7.0, pgCC-5.0
  // IBM xlC has correct large file support
  // Remove this ifdef when other compilers catch up...
  FILE* outfile;
  if ( ctxt_.onpe0() )
  {
    outfile = fopen(filename_cstr,"w");
    string header("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<fpmd:sample xmlns:fpmd=\"");
    header += qbox_xmlns();
    header += string("\"\n");
    header += string(
    " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    header += string(" xsi:schemaLocation=\"");
    header += qbox_xmlns();
    header += string(" sample.xsd\">\n");
    off_t len = header.size();
    fwrite(header.c_str(),sizeof(char),len,outfile);
    
    string desc = string("<description> ") +
      description +
      string(" </description>\n");
    
    ostringstream ss("");
    ss << desc;
    ss << s.atoms;
    string str = ss.str();
    const char* buf = str.c_str();
    len = str.length();
    fwrite(buf,sizeof(char),len,outfile);
  }
  
  if ( !atomsonly )
  {
    s.wf.write(outfile,encoding,"wavefunction");
    if ( s.wfv != 0 )
      s.wfv->write(outfile,encoding,"wavefunction_velocity");
  }
  
  if ( ctxt_.onpe0() )
  {
    char *trailer = "</fpmd:sample>\n";
    fwrite(trailer,sizeof(char),strlen(trailer),outfile);
    fclose(outfile);
  }
#else
  ofstream os;
  if ( ctxt_.onpe0() )
  {
    os.open(filename_cstr);
    cout << "  <!-- SaveCmd: saving to file " << filename 
         << ", encoding=" << encoding << " -->" << endl;
    
    os 
<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
<<"<fpmd:sample xmlns:fpmd=\""
<< qbox_xmlns()
<< "\"\n" 
<<" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
<<" xsi:schemaLocation=\""
<< qbox_xmlns() << " sample.xsd\">"
<< endl;

    os << "<description> " << description     
       << " </description>" << endl;
    os << s.atoms;
  }
    
  if ( !atomsonly )
  {
    s.wf.print(os,encoding,"wavefunction");
    if ( s.wfv != 0 )
      s.wfv->print(os,encoding,"wavefunction_velocity");
  }

  if ( ctxt_.onpe0() )
    os << "</fpmd:sample>" << endl;   

  os.close();
#endif
}
