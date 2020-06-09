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
// testjacobi.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Test the Jacobi implementation of the Matrix class
//
// use: testjacobi nprow npcol n nb
//

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <algorithm>

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Context.h"
#include "Matrix.h"
#include "jacobi.h"
using namespace std;

double aa(int i, int j) { return 1.0/(i+1)+2.0*i/(j+1); }
double frank(int n, int i, int j) { return n - max(i,j); }
//double frank(int n, int i, int j) { return n - fabs(i-j); }
double bb(int i, int j) { return i-j-3; }

int main(int argc, char **argv)
{
  // use: testjacobi nprow npcol n nb
  int mype;
  int npes;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
#else
  npes=1;
  mype=0;
#endif

  //char* infilename = argv[1];
  //ifstream infile(infilename);

  Timer tm;
  if ( argc != 5 )
  {
    cout << "use: testjacobi nprow npcol n nb" << endl;
    return 1;
  }
  int nprow=atoi(argv[1]);
  int npcol=atoi(argv[2]);
  int m_a=atoi(argv[3]);
  int mb_a=atoi(argv[4]);
  int n_a = m_a;
  int nb_a = mb_a;
  if(mype == 0)
  {
    //infile >> nprow >> npcol;
    cout<<"nprow="<<nprow<<", npcol="<<npcol<<endl;
    //infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
  }
#ifdef USE_MPI
  MPI_Bcast(&nprow, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcol, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  {
    Context ctxt(MPI_COMM_WORLD,nprow,npcol);

    if ( mype == 0 )
    {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }

    DoubleMatrix a(ctxt,m_a,n_a,mb_a,nb_a);
    DoubleMatrix u(ctxt,m_a,n_a,mb_a,nb_a);
    vector<double> e(n_a);

    if ( mype == 0 )
    {
      cout << " m_a x n_a / mb_a x nb_a / ta = "
           << a.m() << "x" << a.n() << " / "
           << a.mb() << "x" << a.nb() << endl;
    }

    for ( int m = 0; m < a.nblocks(); m++ )
      for ( int l = 0; l < a.mblocks(); l++ )
        for ( int y = 0; y < a.nbs(m); y++ )
          for ( int x = 0; x < a.mbs(l); x++ )
          {
            int i = a.i(l,x);
            int j = a.j(m,y);
            //double aij = a.i(l,x) * 10 + a.j(m,y);
            double aij = frank(a.n(),i,j);
            int iii = x + l*a.mb();
            int jjj = y + m*a.nb();
            int ival = iii + jjj * a.mloc();
            a[ival] = aij;
          }
    //cout << a;

    const int maxsweep = 30;
    tm.start();
    double tol = 1.e-10;
    int nsweep = jacobi(maxsweep,tol,a,u,e);
    tm.stop();
    if ( ctxt.onpe0() ) cout << " nsweep: " << nsweep << endl;
    if (ctxt.mype() == 0) cout << "Jacobi time: " << tm.real() << endl;

    sort(e.begin(),e.end());

    vector<double> e_exact(e.size());
    for ( int i = 0; i < e_exact.size(); i++ )
    {
      e_exact[i] = 1.0 / ( 2.0 * ( 1.0 - cos( ((2*i+1)*M_PI)/(2*a.n()+1)) ) );
    }
    sort(e_exact.begin(),e_exact.end());

    double asum = 0.0;
    for ( int i = 0; i < e.size(); i++ )
    {
      //cout << ctxt.mype() << ": eig[" << i << "]= " << e[i]
      //     << "  " << e_exact[i] << endl;
      asum += fabs(e[i]-e_exact[i]);
    }
    cout << " abs eigenvalue error: " << asum << endl;
    //cout << " u=" << endl;
    //cout << u;
    //cout << " au=" << endl;
    //cout << a;

#if 1
    // Eigenvalues and eigenvectors of c if c is square
    if ( a.m() == a.n() && a.mb() == a.nb() )
    {
      for ( int m = 0; m < a.nblocks(); m++ )
        for ( int l = 0; l < a.mblocks(); l++ )
          for ( int y = 0; y < a.nbs(m); y++ )
            for ( int x = 0; x < a.mbs(l); x++ )
            {
              int i = a.i(l,x);
              int j = a.j(m,y);
              int iii = x + l*a.mb();
              int jjj = y + m*a.nb();
              int ival = iii + jjj * a.mloc();
              if ( i == j )
                a[ival] = i + 1.e-5*drand48();
              else
                a[ival] = 1.e-5*drand48();
            }
      tm.reset();
      tm.start();
      if (mype == 0) cout << "Eigenproblem... ";
      DoubleMatrix z(a.context(),a.n(),a.n(),a.nb(),a.nb());
      valarray<double> w(a.m());
      a.syev('l',w,z);
      //c.syevx('l',w,z,1.e-5);
      if (mype == 0) cout << " done" << endl;
      tm.stop();
      if (mype == 0) cout << "Eigenproblem time: " << tm.real() << endl;
    }
#endif
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
}
