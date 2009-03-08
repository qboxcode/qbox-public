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
// testAndersonMixer.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testAndersonMixer.C,v 1.7 2009-03-08 01:16:33 fgygi Exp $


#include <iostream>
#include <vector>
using namespace std;

#include "Context.h"
#include "AndersonMixer.h"

// use: testAndersonMixer ndim nmax niter
int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    if ( argc != 4 )
    {
      cout << " use: testAndersonMixer ndim nmax niter" << endl;
      return 1;
    }
    const int ndim = atoi(argv[1]);
    const int nmax = atoi(argv[2]);
    const int niter = atoi(argv[3]);

    Context ctxt;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    // cout << " Process " << ctxt.mype() << " on " << processor_name << endl;

    const double alpha = 0.1;
    vector<double> x,f,xbar,fbar;
    x.resize(ndim);
    f.resize(ndim);
    xbar.resize(ndim);
    fbar.resize(ndim);

    AndersonMixer mixer(ndim,nmax,&ctxt);

    for ( int i = 0; i < ndim; i++ )
      x[i] = (i+5);

    for ( int iter = 0; iter < niter; iter++ )
    {
      // compute gradient
#if 1
      // quadratic form
      for ( int i = 0; i < ndim; i++ )
        f[i] = -(i+1) * (x[i]-i);
        //f[i] = -(0.1*i+1)*(x[i]-i);
#if 1
      // precondition f
      for ( int i = 0; i < ndim; i++ )
        f[i] *= 1.0/(i+3);
#endif

#else
      // Rosenbrock function
      assert(ndim%2==0);
      // f = sum_i^(n-1) (1-x_i)^2 + 100*(x_i+1 - x_i^2)^2
      f[0] = 0.0;
      for ( int i = 0; i < ndim-1; i++ )
      {
        f[i+1] = 200 * ( x[i+1] - x[i]*x[i] );
        f[i] += -2.0 * ( 1.0 - x[i] ) + 200 * x[i] * ( x[i+1] - x[i]*x[i] );
      }
#endif

      double resnorm = 0.0;
      for ( int i = 0; i < ndim; i++ )
        resnorm += f[i]*f[i];
#if USE_MPI
      ctxt.dsum(1,1,&resnorm,1);
#endif

      cout << " resnorm: " << sqrt(resnorm) << endl;
      mixer.update(&x[0],&f[0],&xbar[0],&fbar[0]);

#if 0
      // precondition fbar
      for ( int i = 0; i < ndim; i++ )
        fbar[i] /= (i+3);
        //fbar[i] *= 1.0/(i+2);
#endif

      for ( int i = 0; i < ndim; i++ )
        x[i] = xbar[i] + alpha * fbar[i];
    }
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
