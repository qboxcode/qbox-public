////////////////////////////////////////////////////////////////////////////////
//
// LoadCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoadCmd.C,v 1.9 2008-02-12 05:39:18 fgygi Exp $

#include "LoadCmd.h"
#include "SampleReader.h"
#include "Sample.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int LoadCmd::action(int argc, char **argv)
{
  if ( (argc != 2 && argc !=3) && ui->onpe0() )
  {
    cout << "  use: load [-serial] uri" << endl;
    return 1;
  }

  int iarg = 1;
  bool serial = false;

  if ( !strcmp(argv[iarg],"-serial") )
  {
    serial = true;
    iarg++;
  }

  if ( ui->onpe0() )
    cout << " LoadCmd: loading from " << argv[iarg] << endl;

  SampleReader s_reader(s->ctxt_);

  try
  {
    s_reader.readSample(*s,argv[iarg],serial);
  }
  catch ( const SampleReaderException& e )
  {
    cout << " SampleReaderException caught in LoadCmd:" << endl;
    cout << e.msg << endl;
  }
  catch (...)
  {
    cout << " LoadCmd: cannot load Sample" << endl;
  }

  s->ctxt_.barrier();

  return 0;
}
