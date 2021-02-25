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
// testSampleReader.cpp
//
// Test functionality of SampleReader
// use: ./testSampleReader sample.xml
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
using namespace std;

#include "MPIdata.h"
#include "Context.h"
#include "SlaterDet.h"
#include "UnitCell.h"
#include "Sample.h"
#include "D3vector.h"
#include "SampleReader.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int ngb = size;
  int nstb = 1;
  int nkpb = 1;
  int nspb = 1;
  MPIdata::set(ngb,nstb,nkpb,nspb);
  {
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << MPIdata::rank() << " on " << processor_name << endl;

    Sample* s = new Sample;

    SampleReader s_reader;
    bool serial = false;
    if ( argc != 2 )
    {
      cout << "use: testSampleReader {file|URI}" << endl;
      return 1;
    }
    const char* filename = argv[1];

    try
    {
      s_reader.readSample(*s,filename,serial);
    }
    catch ( const SampleReaderException& e )
    {
      cout << " SampleReaderException caught:" << endl;
      cout << e.msg << endl;
    }
    catch (...)
    {
      cout << " testSampleReader: cannot load Sample" << endl;
    }

    if ( MPIdata::onpe0() )
    {
      cout << filename << endl;
      cout << s->atoms.size() << " atoms" << endl;
      cout << s->wf.nel() << " electrons" << endl;
    }
  }
  MPI_Finalize();
  return 0;
}
