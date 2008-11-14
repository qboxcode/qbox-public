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
// $Id: testAndersonMixer.C,v 1.6 2008-11-14 04:08:05 fgygi Exp $


#include <iostream>
#include <vector>
using namespace std;

#include "Context.h"
#include "AndersonMixer.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    const int ndim = 20;
    const int niter = 10;

    Context ctxt;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    // cout << " Process " << ctxt.mype() << " on " << processor_name << endl;

    const double alpha = 1.0;
    vector<double> x,f,xlast,fbar;
    double theta;
    x.resize(ndim);
    xlast.resize(ndim);
    f.resize(ndim);
    fbar.resize(ndim);

    AndersonMixer mixer(x.size(),&ctxt);

    for ( int i = 0; i < ndim; i++ )
      x[i] = (i+5);

    xlast = x;

    for ( int iter = 0; iter < niter; iter++ )
    {
      // compute gradient
#if 1
      // quadratic form
      for ( int i = 0; i < ndim; i++ )
        //f[i] = -(i+1) * (x[i]-1.0) - 0.1 * x[i]*x[i]*x[i];
        f[i] = -(i+1) * (x[i]-i);
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
      for ( int i = 0; i < x.size(); i++ )
        resnorm += f[i]*f[i];
#if USE_MPI
      ctxt.dsum(1,1,&resnorm,1);
#endif


      mixer.update(&f[0],&theta,&fbar[0]);
      cout << " theta: " << theta << " resnorm: " << resnorm << endl;

      for ( int i = 0; i < x.size(); i++ )
      {
        // xbar = x + theta * ( x - xlast )
        // x = xbar + fbar
        double xtmp = x[i];
        x[i] += theta * ( x[i] - xlast[i] ) + alpha * fbar[i];
        xlast[i] = xtmp;
      }
    }
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
