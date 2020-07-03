////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2020 The Regents of the University of California
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
// testMPIdata.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Test of the MPIdata class
//
// use: testMPIdata ngb nstb nspb nkpb

#include "MPIdata.h"
#include<iostream>
using namespace std;

int main(int argc, char** argv)
{
  if ( argc != 5 )
  {
    cerr << "use: testMPIdata ngb nstb nspb nkpb" << endl;
    return 1;
  }
  MPI_Init(&argc,&argv);
  int ngb = atoi(argv[1]);
  int nstb = atoi(argv[2]);
  int nspb = atoi(argv[3]);
  int nkpb = atoi(argv[4]);

  MPIdata::set(ngb,nstb,nspb,nkpb);

  int npes;
  MPI_Comm_size(MPIdata::comm(),&npes);
  cout << " rank=" << MPIdata::rank() << "    comm size=" << npes << endl;
  MPI_Comm_size(MPIdata::g_comm(),&npes);
  cout << " rank=" << MPIdata::rank() << "  g_comm size=" << npes << endl;
  MPI_Comm_size(MPIdata::st_comm(),&npes);
  cout << " rank=" << MPIdata::rank() << " st_comm size=" << npes << endl;
  MPI_Comm_size(MPIdata::sp_comm(),&npes);
  cout << " rank=" << MPIdata::rank() << " sp_comm size=" << npes << endl;
  MPI_Comm_size(MPIdata::kp_comm(),&npes);
  cout << " rank=" << MPIdata::rank() << " kp_comm size=" << npes << endl;
  MPI_Comm_size(MPIdata::sd_comm(),&npes);
  cout << " rank=" << MPIdata::rank() << " sd_comm size=" << npes << endl;

  MPI_Finalize();
}
