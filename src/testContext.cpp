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
// testContext.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <vector>
#include <cstdlib>
using namespace std;
#include "Context.h"

int main(int argc, char **argv)
{
  int mype;
  int size;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  if ( argc != 5 )
  {
    cout << "use: testContext ngb nstb nspb nkpb" << endl;
    return 1;
  }
  // read number of G-vector blocks, states blocks, spin blocks, kp blocks
  int ngb = atoi(argv[1]);
  int nstb = atoi(argv[2]);
  int nspb = atoi(argv[3]);
  int nkpb = atoi(argv[4]);
  // check that all numbers of blocks are positive
  assert(ngb>0);
  assert(nstb>0);
  assert(nspb>0);
  assert(nkpb>0);

  // create 4-dim Cart_comm
  assert(ngb*nstb*nspb*nkpb == size);

  int ndims=4;
  int dims[4] = { ngb, nstb, nspb, nkpb };
  // plan for cyclic rotation of states blocks
  int periods[4] = { 0, 1, 0, 0};
  int reorder = 0;

  MPI_Comm comm;
  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm);

  // print 4-dim coordinates of each task
  int coords[4];
  MPI_Cart_coords(comm,mype,4,coords);

  cout << mype << ": coords: "
       << coords[0] << " " << coords[1] << " "
       << coords[2] << " " << coords[3] << endl;

  { // start Context scope

    Context ctxt(comm);
    for ( int i = 0; i < size; i++ )
    {
      MPI_Barrier(comm);
      if ( i == mype )
        cout << mype << ":" << ctxt.mype() << ":" << ctxt.myproc()
         << " base: " << ctxt;
    }
    vector<Context*> c;

    c.push_back(new Context(comm,ngb,nstb));
    cout << ctxt.mype() << ": " << *c[0];

    if ( c[0]->active() )
      cout << ctxt.mype() << ": c[0] is active" << endl;

    // test dgsum2d function
    // add along rows, then along columns
    double a = c[0]->mype();
    cout << c[0]->mype() << ": a     = " << a << endl;
    c[0]->dsum('R',1,1,&a,1);
    cout << c[0]->mype() << ": a_sum_row = " << a << endl;
    c[0]->dsum('C',1,1,&a,1);
    cout << c[0]->mype() << ": a_sum_all = " << a << endl;

    for ( int i = 0; i < c.size(); i++ )
    {
      delete c[i];
    }

  } // end Context scope

  MPI_Comm_free(&comm);
  MPI_Finalize();
}
