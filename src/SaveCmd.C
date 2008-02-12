////////////////////////////////////////////////////////////////////////////////
//
// SaveCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SaveCmd.C,v 1.14 2008-02-12 05:39:18 fgygi Exp $


#include "SaveCmd.h"
#include "SampleWriter.h"
#include "isodate.h"
#include "release.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SaveCmd::action(int argc, char **argv)
{
  string usage("  Use: save [-text|-base64] [-atomsonly] [-serial] filename");
  if ( !(argc>=2 && argc<=4 ) )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }

  // set default encoding
  bool base64 = true;
  bool atomsonly = false;
  bool serial = false;
  char* filename = 0;

  // check for -text or -base64 or -atomsonly or -serial arguments
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
    else if ( arg=="-serial" )
    {
      serial = true;
    }
    else if ( arg[0] != '-' && i == argc-1 )
    {
      filename = argv[i];
    }
    else
    {
      if ( ui->onpe0() )
        cout << usage << endl;
      return 1;
    }
  }

  if ( filename == 0 )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }
  SampleWriter swriter(s->ctxt_);
  string description = string(" Created ") + isodate() +
                       string(" by qbox-") + release() + string(" ");
  swriter.writeSample(*s, filename, description, base64, atomsonly, serial);

  return 0;
}
