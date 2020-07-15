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
// use: testMPIdata ngb nstb nkpb nspb

#include "MPIdata.h"
#include<iostream>
#include<cstdlib> // atoi
using namespace std;

int main(int argc, char** argv)
{
  if ( argc != 5 )
  {
    cerr << "use: testMPIdata ngb nstb nkpb nspb" << endl;
    return 1;
  }
  MPI_Init(&argc,&argv);
  int ngb = atoi(argv[1]);
  int nstb = atoi(argv[2]);
  int nkpb = atoi(argv[3]);
  int nspb = atoi(argv[4]);

  MPIdata::set(ngb,nstb,nkpb,nspb);

  cout << " rank=" << MPIdata::rank() << " igb / ngb = "
       << MPIdata::igb() << " / " << MPIdata::ngb() << endl;
  cout << " rank=" << MPIdata::rank() << " istb / nstb = "
       << MPIdata::istb() << " / " << MPIdata::nstb() << endl;
  cout << " rank=" << MPIdata::rank() << " ikpb / nkpb = "
       << MPIdata::ikpb() << " / " << MPIdata::nkpb() << endl;
  cout << " rank=" << MPIdata::rank() << " ispb / nspb = "
       << MPIdata::ispb() << " / " << MPIdata::nspb() << endl;

  int coords[2];
  MPI_Cart_coords(MPIdata::sd_comm(),MPIdata::sd_rank(),2,coords);
  cout << " rank=" << MPIdata::rank() << " sd_rank = "
       <<MPIdata::sd_rank()
       << " coords=(" << coords[0] << "," << coords[1] << ")" << endl;

  MPI_Barrier(MPIdata::sd_comm());

  MPI_Finalize();
}
