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
// testBasis.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Basis.h"
#include <iostream>
#include <new>
#include <cstdlib>
#include <cassert>
using namespace std;

#include <mpi.h>

int main(int argc, char **argv)
{
  // use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z ecut kx ky kz npr npc
  MPI_Init(&argc,&argv);
  {
    if ( argc !=16 )
    {
      cout << " use: testBasis a0x a0y a0z a1x a1y a1z a2x a2y a2z"
           << " ecut kx ky kz npr npc" << endl;
      return 1;
    }
    const D3vector a0(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    const D3vector a1(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    const D3vector a2(atof(argv[7]),atof(argv[8]),atof(argv[9]));

    double ecut = atof(argv[10]);
    D3vector kpoint(atof(argv[11]),atof(argv[12]),atof(argv[13]));
    int npr = atoi(argv[14]);
    int npc = atoi(argv[15]);

    UnitCell cell(a0,a1,a2);

    // create cartesian communicator
    int ndims=2;
    int dims[2] = {npr, npc};
    int periods[2] = { 0, 0};
    int reorder = 0;
    MPI_Comm cartcomm;
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&cartcomm);

    // partition the cartesian communicator in columns
    MPI_Comm colcomm;
    int remain_dims[2] = { 1, 0 };
    MPI_Cart_sub(cartcomm,remain_dims,&colcomm);

    Basis basis(colcomm,kpoint);
    try
    {
      basis.resize(cell,cell,ecut);
    }
    catch ( bad_alloc )
    {
      cout << " bad_alloc caught in Basis::resize" << endl;
      throw;
    }

    int npes, mype;
    MPI_Comm_size(colcomm,&npes);
    MPI_Comm_rank(colcomm,&mype);
    for (int ipe = 0; ipe < npes; ipe++ )
    {
      MPI_Barrier(colcomm);
      if ( ipe == mype )
      {
        cout << basis;
        cout << endl;
      }
    }

    //Basis b2(basis);
    //cout << b2;
  }
  MPI_Finalize();
}
