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
// testContext.C
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
  int npes;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  int nr = atoi(argv[1]);
  int nc = atoi(argv[2]);

  { // start Context scope

    Context ctxt(MPI_COMM_WORLD);
    for ( int i = 0; i < npes; i++ )
    {
      MPI_Barrier(MPI_COMM_WORLD);
      if ( i == mype )
        cout << mype << ":" << ctxt.mype() << ":" << ctxt.myproc()
         << " base: " << ctxt;
    }
    vector<Context*> c;

    c.push_back(new Context(MPI_COMM_WORLD,nr,nc));
    cout << ctxt.mype() << ": " << *c[0];

    if ( *c[0] )
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

  MPI_Finalize();
}
