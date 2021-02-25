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
// testSample.cpp
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

    Sample s;

    D3vector a(18, 0, 0);
    D3vector b( 0,18, 0);
    D3vector c( 0, 0,18);
    UnitCell uc(a,b,c);
    double ecut = 25.0;
    s.wf.resize(uc,uc,ecut);
    s.wf.set_nel(12*54);

    s.wf.randomize(1.e-4);
    s.wf.gram();
    cout << " ortho_error: " << s.wf.sd(0,0)->ortho_error() << endl;
  }
  MPI_Finalize();
  return 0;
}
