////////////////////////////////////////////////////////////////////////////////
//
// LoadCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: LoadCmd.C,v 1.6 2004-09-14 22:24:11 fgygi Exp $

#include "LoadCmd.h"
#include "SampleReader.h"
#include "Sample.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int LoadCmd::action(int argc, char **argv)
{
  if ( (argc != 2 && argc !=3) && ui->onpe0() )
  {
    cout << "  <!-- use: load [-serial] uri -->" << endl;
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
    cout << " <!-- LoadCmd: loading from " << argv[iarg] << " -->" << endl;

  SampleReader s_reader(s->ctxt_);
  
  if ( ui->onpe0() )
    cout << " <!--" << endl;
    
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
  
  // If only <atomset> was read, set nel for the wavefunction
  //cout << " LoadCmd: atoms.nel() = " << s->atoms.nel() << endl;
  //cout << " LoadCmd: wf.nel() =    " << s->wf.nel() << endl;
  if ( s->wf.nel() != s->atoms.nel() )
  {
    s->wf.set_nel(s->atoms.nel());
    s->wf.update_occ(0.0);
  }
  //cout << " LoadCmd: atoms.nel() = " << s->atoms.nel() << endl;
  //cout << " LoadCmd: wf.nel() =    " << s->wf.nel() << endl;

  if ( ui->onpe0() )
    cout << " -->" << endl;
    
  return 0;
}
