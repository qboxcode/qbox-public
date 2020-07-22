////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
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
// LoadCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

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

  // Reset current sample
  // cout << "atomset before reset: nsp: " << s->atoms.nsp() << endl;
  s->reset();
  // cout << "atomset after reset: nsp: " << s->atoms.nsp() << endl;

  SampleReader s_reader;

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

  MPI_Barrier(MPIdata::comm());

  return 0;
}
