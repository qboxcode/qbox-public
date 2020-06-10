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
// testjade.cpp
//
// use: testjade nprow npcol m mb
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <algorithm>
using namespace std;

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Context.h"
#include "Matrix.h"
#include "jade.h"

const bool print = false;

enum MatrixType { FRANK, SCALED_FRANK, LAPLACIAN, SMALL };
const MatrixType mtype = SCALED_FRANK;

double aa(int i, int j) { return 1.0/(i+1)+2.0*i/(j+1); }
double bb(int i, int j) { return i-j-3; }

double frank(int n, int i, int j) { return n - max(i,j); }
double aijf(int n, int k, int i, int j)
{
  switch ( mtype )
  {
    case FRANK :
      if ( k == 0 )
        return frank(n,i,j);
      else
        return frank(n,n-i,n-j);

    case SCALED_FRANK :
      if ( k == 0 )
        return frank(n,i,j)/(n*n);
      else
        return frank(n,n-i,n-j)/(n*n);

    case LAPLACIAN :
      if ( i == j )
      {
        return 2.0 + k;
      }
      else if ( (i-j)*(i-j) == 1 )
      {
        return 1.0;
      }

    case SMALL :
    {
      if ( k == 0 )
      {
        if ( i == j && i == 0 )
          return 1.0;
        else if ( i == j && i > 0 )
          return 4.0;
        else
          return 0.0;
      }
      else
      {
        if ( i == j && i == 0 )
          return 1.0;
        else if ( i == j && i > 0 )
          return 1.0;
        else
          return 0.5;
      }
    }
  }
  return 0.0;
}

int main(int argc, char **argv)
{
  // use: testjade nprow npcol n nb
  const int maxsweep = 30;
  const double tol = 1.e-8;
  const int nmat = 2;

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

  Timer tm;
  if ( argc != 5 )
  {
    cout << "use: testjade nprow npcol n nb" << endl;
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
    cout << "nprow=" << nprow << ", npcol=" << npcol << endl;
    //infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout << "nmat=" << nmat << " m_a=" << m_a << ", n_a=" << n_a << endl;
    cout << "tol=" << tol << endl;
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

    vector<DoubleMatrix*> a(nmat);
    vector<vector<double> > adiag(nmat);
    for ( int k = 0; k < nmat; k++ )
    {
      a[k] = new DoubleMatrix(ctxt,m_a,n_a,mb_a,nb_a);
      adiag[k].resize(n_a);
    }
    DoubleMatrix u(ctxt,m_a,n_a,mb_a,nb_a);

    if ( mype == 0 )
    {
      cout << " m_a x n_a / mb_a x nb_a / ta = "
           << a[0]->m() << "x" << a[0]->n() << " / "
           << a[0]->mb() << "x" << a[0]->nb() << endl;
    }

    for ( int k = 0; k < nmat; k++ )
    {
      for ( int m = 0; m < a[k]->nblocks(); m++ )
        for ( int l = 0; l < a[k]->mblocks(); l++ )
          for ( int y = 0; y < a[k]->nbs(m); y++ )
            for ( int x = 0; x < a[k]->mbs(l); x++ )
            {
              int i = a[k]->i(l,x);
              int j = a[k]->j(m,y);
              //double aij = a.i(l,x) * 10 + a.j(m,y);
              double aij = aijf(a[k]->n(),k,i,j);
              int iii = x + l*a[k]->mb();
              int jjj = y + m*a[k]->nb();
              int ival = iii + jjj * a[k]->mloc();
              (*a[k])[ival] = aij;
            }
      //cout << " a[" << k << "]=" << endl;
      //cout << (*a[k]);
    }

    tm.start();
    int nsweep = jade(maxsweep,tol,a,u,adiag);
    tm.stop();
    if ( ctxt.onpe0() )
    {
      cout << " m=" << m_a << " mb=" << mb_a
           << " ctxt: " << ctxt.nprow() << "x" << ctxt.npcol()
           << " nsweep=" << nsweep << " time: " << tm.real() << endl;
    }

    for ( int k = 0; k < nmat; k++ )
      sort(adiag[k].begin(),adiag[k].end());

    if ( nmat == 1 && (mtype == FRANK || mtype == SCALED_FRANK) )
    {
      vector<double> e_exact(adiag[0].size());
      for ( int i = 0; i < e_exact.size(); i++ )
      {
        e_exact[i] = 1.0 / ( 2.0 * ( 1.0 - cos( ((2*i+1)*M_PI)/(2*n_a+1)) ) );
        if ( mtype == SCALED_FRANK )
        {
          int n = a[0]->n();
          e_exact[i] /= (n*n);
        }
      }
      sort(e_exact.begin(),e_exact.end());

      if ( mype == 0 )
      {
        for ( int k = 0; k < nmat; k++ )
        {
          double asum = 0.0;
          for ( int i = 0; i < e_exact.size(); i++ )
          {
            //cout << ctxt.mype() << ": eig[" << k << "][" << i << "]= "
            //     << adiag[k][i]
            //     << "  " << e_exact[i] << endl;
            asum += fabs(adiag[k][i]-e_exact[i]);
          }
          cout << "a[" << k << "] sum of abs eigenvalue errors: "
               << asum << endl;
        }
      }

    }

    if ( print )
    {
      // a[k] contains AU.
      // compute the product C = U^T A U
      DoubleMatrix c(*a[0]);
      for ( int k = 0; k < nmat; k++ )
      {
        c.gemm('t','n',1.0,u,(*a[k]),0.0);
        if ( mype == 0 )
        {
          cout << " a[" << k << "]=" << endl;
          cout << c;
        }
      }
      cout << " u=" << endl;
      cout << u;
    }
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
}
