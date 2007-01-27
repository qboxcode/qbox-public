////////////////////////////////////////////////////////////////////////////////
//
// SaveCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveCmd.C,v 1.11 2007-01-27 23:49:41 fgygi Exp $


#include "SaveCmd.h"
#include "SampleWriter.h"
#include "isodate.h"
#include "release.h"

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
  bool base64 = true;
  bool atomsonly = false;
  char* filename = 0;
  
  // check for -text or -base64 or -atomsonly arguments
  for ( int i = 1; i < argc; i++ )
  {
    string arg(argv[i]);
    
    if ( arg=="-text" )
    {
      base64 = false;
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
  SampleWriter swriter(s->ctxt_);
  string description = string(" Created ") + isodate() +
                       string(" by qbox-") + release() + string(" ");
  swriter.writeSample(*s, filename, description, base64, atomsonly);

  return 0;
}
