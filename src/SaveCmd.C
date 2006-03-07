////////////////////////////////////////////////////////////////////////////////
//
// SaveCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveCmd.C,v 1.9 2006-03-07 07:36:57 fgygi Exp $


#include "SaveCmd.h"
#include "fstream"
#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"

#ifdef USE_CSTDIO_LFS
#include <cstdio>
#include <cstdlib>
#include <sstream>
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SaveCmd::action(int argc, char **argv)
{
  if ( !(argc>=2 && argc<=4 ) )
  {
    if ( ui->onpe0() )
      cout << "  <!-- use: save [-text|-base64] [-atomsonly] filename -->" 
           << endl;
    return 1;
  }
  
  // set default encoding
  string encoding = "base64";
  bool atomsonly = false;
  char* filename = 0;
  
  // check for -text or -base64 or -atomsonly arguments
  for ( int i = 1; i < argc; i++ )
  {
    string arg(argv[i]);
    
    if ( arg=="-text" )
    {
      encoding = "text";
    }
    else if ( arg=="-base64" )
    {
      encoding = "base64";
    }
    else if ( arg=="-atomsonly" )
    {
      atomsonly = true;
    }
    else if ( arg[0] != '-' && i == argc-1 )
    {
      filename = argv[i];
    }
    else
    {
      if ( ui->onpe0() )
        cout << "  <!-- use: save [-text|-base64] [-atomsonly] filename -->" 
             << endl;
      return 1;
    }
  }
  
  if ( filename == 0 )
  {
    if ( ui->onpe0() )
      cout << "  <!-- use: save [-text|-base64] [-atomsonly] filename -->" 
           << endl;
    return 1;
  }

#if USE_CSTDIO_LFS
  // This section for compilers with broken large file support
  // As of 2003-09-30, this includes gcc-3.2, icc-7.0, pgCC-5.0
  // IBM xlC has correct large file support
  // Remove this ifdef when other compilers catch up...
  FILE* outfile;
  if ( ui->onpe0() )
  {
    outfile = fopen(filename,"w");
    string header("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<fpmd:sample xmlns:fpmd=\"");
    header += qbox_xmlns();
    header += string("\"\n");
    header += string(
    " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
    " xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
    header += string(" xsi:schemaLocation=\"");
    header += qbox_xmlns();
    header += string(" sample.xsd\">\n");
    off_t len = header.size();
    fwrite(header.c_str(),sizeof(char),len,outfile);
    
    string desc = string("<description> Created ") +
      isodate() +string(" by qbox-") + release() +
      string(" </description>\n");
    
    ostringstream ss("");
    ss << desc;
    ss << s->atoms;
    string str = ss.str();
    const char* buf = str.c_str();
    len = str.length();
    fwrite(buf,sizeof(char),len,outfile);
  }
  
  if ( !atomsonly )
  {
    s->wf.write(outfile,encoding,"wavefunction");
    if ( s->wfv != 0 )
      s->wfv->write(outfile,encoding,"wavefunction_velocity");
  }
  
  if ( ui->onpe0() )
  {
    char *trailer = "</fpmd:sample>\n";
    fwrite(trailer,sizeof(char),strlen(trailer),outfile);
    fclose(outfile);
  }
#else
  ofstream os;
  if ( ui->onpe0() )
  {
    os.open(filename);
    cout << "  <!-- SaveCmd: saving to file " << filename 
         << ", encoding=" << encoding << " -->" << endl;
    
    os 
<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
<<"<fpmd:sample xmlns:fpmd=\""
<< qbox_xmlns()
<< "\"\n" 
<<" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
<<" xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
<<" xsi:schemaLocation=\""
<< qbox_xmlns() << " sample.xsd\">"
<< endl;

    os << "<description> Created " << isodate() << " by qbox-" << release()
       << " </description>" << endl;
    os << s->atoms;
  }
    
  if ( !atomsonly )
  {
    s->wf.print(os,encoding,"wavefunction");
    if ( s->wfv != 0 )
      s->wfv->print(os,encoding,"wavefunction_velocity");
  }

  if ( ui->onpe0() )
    os << "</fpmd:sample>" << endl;   

  os.close();
#endif
  
  return 0;
}
