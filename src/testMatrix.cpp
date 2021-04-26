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
// testMatrix.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testMatrix input_file [-check]
// input_file:
// nprow npcol
// m_a n_a mb_a nb_a transa
// m_b n_b mb_b nb_b transb
// m_c n_c mb_c nb_c
//

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
using namespace std;

#include "MPIdata.h"
#include "Timer.h"

#include <mpi.h>
#include "Context.h"
#include "Matrix.h"

double aa(int i, int j) { return 1.0/(i+1)+2.0*i/(j+1); }
double bb(int i, int j) { return i-j-3; }

int main(int argc, char **argv)
{
  int mype;
  int npes;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  char* infilename = argv[1];
  ifstream infile(infilename);

  assert(argc == 2 || argc == 3);
  if ( !(argc==2 || argc==3) )
  {
    cout << "use: testMatrix inputfile [-check]" << endl;
    return 1;
  }
  bool tcheck = false;
  if ( argc == 3 )
  {
    if ( !strcmp(argv[2],"-check") )
      tcheck = true;
    else
    {
      cerr << " invalid argv[2]" << endl;
      MPI_Abort(MPI_COMM_WORLD,2);
    }
  }
  Timer tm;
  int nprow, npcol;
  int m_a, n_a, mb_a, nb_a;
  int m_b, n_b, mb_b, nb_b;
  int m_c, n_c, mb_c, nb_c;
  char ta, tb;
  if(mype == 0)
  {
    infile >> nprow >> npcol;
    cout<<"nprow="<<nprow<<", npcol="<<npcol<<endl;
    infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
    infile >> m_b >> n_b >> mb_b >> nb_b >> tb;
    cout<<"m_b="<<m_b<<", n_b="<<n_b<<endl;
    infile >> m_c >> n_c >> mb_c >> nb_c;
    cout<<"m_c="<<m_c<<", n_c="<<n_c<<endl;
  }
  MPI_Bcast(&nprow, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcol, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ta, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tb, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

  MPIdata::set(nprow,npcol,1,1);

  {
    if ( ta == 'N' ) ta = 'n';
    if ( tb == 'N' ) tb = 'n';

    Context ctxt(MPI_COMM_WORLD,nprow,npcol);

    if ( mype == 0 )
    {
      cout << " Context " << ctxt.ictxt()
           << ": " << ctxt.nprow() << "x" << ctxt.npcol() << endl;
    }

    DoubleMatrix a(ctxt,m_a,n_a,mb_a,nb_a);
    DoubleMatrix b(ctxt,m_b,n_b,mb_b,nb_b);
    DoubleMatrix c(ctxt,m_c,n_c,mb_c,nb_c);

    if ( mype == 0 )
    {
      cout << " m_a x n_a / mb_a x nb_a / ta = "
           << a.m() << "x" << a.n() << " / "
           << a.mb() << "x" << a.nb() << " / " << ta << endl;
      cout << " m_b x n_b / mb_b x nb_b / tb = "
           << b.m() << "x" << b.n() << " / "
           << b.mb() << "x" << b.nb() << " / " << tb << endl;
      cout << " m_c x n_c / mb_c x nb_c      = "
           << c.m() << "x" << c.n() << " / "
           << c.mb() << "x" << c.nb() << endl;
    }

    for ( int m = 0; m < a.nblocks(); m++ )
      for ( int l = 0; l < a.mblocks(); l++ )
        for ( int y = 0; y < a.nbs(m); y++ )
          for ( int x = 0; x < a.mbs(l); x++ )
          {
            int i = a.i(l,x);
            int j = a.j(m,y);
            // double aij = a.i(l,x) * 10 + a.j(m,y);
            double aij = aa(i,j);
            int iii = x + l*a.mb();
            int jjj = y + m*a.nb();
            int ival = iii + jjj * a.mloc();
            a[ival] = aij;
          }

    for ( int m = 0; m < b.nblocks(); m++ )
      for ( int l = 0; l < b.mblocks(); l++ )
        for ( int y = 0; y < b.nbs(m); y++ )
          for ( int x = 0; x < b.mbs(l); x++ )
          {
            int i = b.i(l,x);
            int j = b.j(m,y);
            // double bij = b.i(l,x) * 10 + b.j(m,y);
            double bij = bb(i,j);
            int iii = x + l*b.mb();
            int jjj = y + m*b.nb();
            int ival = iii + jjj * b.mloc();
            b[ival] = bij;
          }

    tm.start();
    double gemm_time = MPI_Wtime();
    c.gemm(ta,tb,1.0,a,b,0.0);
    gemm_time = MPI_Wtime() - gemm_time;
    if ( ctxt.onpe0() )
      cout << "gemm_time: " << gemm_time << endl;
    tm.stop();

    Context csq(MPI_COMM_WORLD,1,1);
    DoubleMatrix ssq(csq,c.n(),c.n(),c.nb(),c.nb());
    double getsub_time = MPI_Wtime();
    ssq.getsub(c,c.m(),c.n(),0,0);
    getsub_time = MPI_Wtime() - getsub_time;
    if ( ctxt.onpe0() )
      cout << "getsub_time: " << getsub_time << endl;

    if ( tcheck )
    {
    cout << " checking results..." << endl;
    for ( int m = 0; m < c.nblocks(); m++ )
      for ( int l = 0; l < c.mblocks(); l++ )
        for ( int y = 0; y < c.nbs(m); y++ )
          for ( int x = 0; x < c.mbs(l); x++ )
          {
            int i = c.i(l,x);
            int j = c.j(m,y);
            double sum = 0.0;
            int kmax = ( ta == 'n' ) ? a.n() : a.m();

            if ( ( ta == 'n' ) && ( tb == 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(i,k) * bb(k,j);
            }
            else if ( ( ta != 'n' ) && ( tb == 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(k,i) * bb(k,j);
            }
            else if ( ( ta == 'n' ) && ( tb != 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(i,k) * bb(j,k);
            }
            else if ( ( ta != 'n' ) && ( tb != 'n' ) )
            {
              for ( int k = 0; k < kmax; k++ )
                sum += aa(k,i) * bb(j,k);
            }

            int iii = x + l*c.mb();
            int jjj = y + m*c.nb();
            int ival = iii + jjj * c.mloc();
            if ( fabs( c[ival] - sum ) > 1.e-8 )
            {
              cout << " error at element (" << i << "," << j << ") "
                   << c[ival] << " " << sum << endl;
              exit(1);
            }
          }
       cout << " results checked" << endl;
    }

    // print dgemm timing
    const int kmax = ( ta == 'n' ) ? a.n() : a.m();
    double buf[2];
    buf[0] = tm.cpu();
    buf[1] = tm.real();
    double cpu, real;
    for ( int irow = 0; irow < nprow; irow++ )
    {
      for ( int icol = 0; icol < npcol; icol++ )
      {
        ctxt.barrier();
        if ( irow != 0 && icol != 0 )
        {
          if ( ctxt.onpe0() )
            ctxt.drecv(1,1,buf,2,irow,icol);
          else
            if ( irow == ctxt.myrow() && icol == ctxt.mycol() )
              ctxt.dsend(1,1,buf,2,0,0);
        }

        if ( ctxt.onpe0() )
        {
          double tcpu = buf[0];
          double treal = buf[1];
          cout << "(" << setw(3) << irow << "," << setw(3) << icol << ") "
               << "CPU/Real: " << setw(8) << tcpu << " / " << setw(8) << treal;
          if ( treal > 0.0 )
            cout << "  MFlops: " << (2.0e-6*m_c*n_c*kmax) / treal;
          cout << endl;
        }
      }
    }

    double norma=a.nrm2();
    if(mype == 0)cout<<"Norm(a)="<<norma<<endl;

    double norm;
    if(mype == 0)cout<<"DoubleMatrix::matgather..."<<endl;
    double*  aa=new double[a.m()*a.n()];
    a.matgather(aa, a.m());
    if(mype == 0)cout<<"DoubleMatrix::init..."<<endl;
    b.init(aa, a.m());
    norm=b.nrm2();
    if ( mype == 0 ) cout << "Norm(b)=" << norm << endl;
    if ( fabs(norm-norma)>0.000001 )
       cout << "DoubleMatrix: problem with matgather/init" << endl;

    if ( c.n() == b.m() && c.m() == b.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::transpose..."<<endl;
      tm.reset();
      tm.start();
      c.transpose(1.0,b,0.0);
      tm.stop();
      if ( mype == 0 ) cout << " transpose time: " << tm.real() << endl;
      norm=c.nrm2();
      if(mype == 0)cout<<"Norm(c)="<<norm<<endl;
    }

    if(mype == 0)cout<<"DoubleMatrix::scal..."<<endl;
    c.scal(0.5);

    if ( a.m() == b.m() && a.n() == b.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::axpy..."<<endl;
      a.axpy(-2., b);
    }

    if ( a.m()==c.m() && a.n()==c.n() && a.mb()==c.mb() && a.nb()==c.nb() )
    {
      if(mype == 0)cout<<"DoubleMatrix::operator=..."<<endl;
      c=a;
    }

    if(mype == 0)cout<<"DoubleMatrix::nrm2..."<<endl;
    norm=c.nrm2();
    if (mype == 0) cout<<"Norm="<<norm<<endl;

    a.identity();
    DoubleMatrix a2(a);
    a -= a2;
    norm = a.nrm2();
    if (mype == 0) cout << "Norm(a)=" << norm << endl;

    // Inverse of c if c is square
    if ( c.m() == c.n() && c.mb() == c.nb() )
    {
      for ( int m = 0; m < c.nblocks(); m++ )
        for ( int l = 0; l < c.mblocks(); l++ )
          for ( int y = 0; y < c.nbs(m); y++ )
            for ( int x = 0; x < c.mbs(l); x++ )
            {
              int i = c.i(l,x);
              int j = c.j(m,y);
              int iii = x + l*c.mb();
              int jjj = y + m*c.nb();
              int ival = iii + jjj * c.mloc();
              if ( i == j )
                c[ival] = i + 1.e-5*drand48();
              else
                c[ival] = 1.e-5*drand48();
            }
      tm.reset();
      tm.start();
      if (mype == 0) cout << "Inverse ... ";
      c.inverse();
      if (mype == 0) cout << " done" << endl;
      tm.stop();
      if (mype == 0) cout << "Inverse time: " << tm.real() << endl;
    }

    // Eigenvalues and eigenvectors of c if c is square
    if ( c.m() == c.n() && c.mb() == c.nb() )
    {
      for ( int m = 0; m < c.nblocks(); m++ )
        for ( int l = 0; l < c.mblocks(); l++ )
          for ( int y = 0; y < c.nbs(m); y++ )
            for ( int x = 0; x < c.mbs(l); x++ )
            {
              int i = c.i(l,x);
              int j = c.j(m,y);
              int iii = x + l*c.mb();
              int jjj = y + m*c.nb();
              int ival = iii + jjj * c.mloc();
              if ( i == j )
                c[ival] = i + 1.e-5*drand48();
              else
                c[ival] = 1.e-5*drand48();
            }
      tm.reset();
      tm.start();
      if (mype == 0) cout << "Eigenproblem... ";
      DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
      valarray<double> w(c.m());
      c.syev('l',w,z);
      //c.syevd('l',w,z);
      //c.syevx('l',w,z,1.e-5);
      if (mype == 0) cout << " done" << endl;
      tm.stop();
      if (mype == 0) cout << "Eigenproblem time: " << tm.real() << endl;

      // complex eigenvalue problem
      ComplexMatrix cc(c.context(),m_c,n_c,mb_c,nb_c);
      for ( int m = 0; m < cc.nblocks(); m++ )
        for ( int l = 0; l < cc.mblocks(); l++ )
          for ( int y = 0; y < cc.nbs(m); y++ )
            for ( int x = 0; x < cc.mbs(l); x++ )
            {
              int i = cc.i(l,x);
              int j = cc.j(m,y);
              int iii = x + l*cc.mb();
              int jjj = y + m*cc.nb();
              int ival = iii + jjj * cc.mloc();
              if ( i == j )
                cc[ival] = i + 1.e-5*drand48();
              else
                cc[ival] = complex<double>(1.e-5*drand48(), 1.e-5*drand48());
            }
      tm.reset();
      tm.start();
      if (mype == 0) cout << "Complex Eigenproblem... ";
      ComplexMatrix zz(cc.context(),cc.n(),cc.n(),cc.nb(),cc.nb());
      valarray<double> ww(cc.m());
      cc.heev('l',ww,zz);
      //c.syevx('l',w,z,1.e-5);
      if (mype == 0) cout << " done" << endl;
      tm.stop();
      if (mype == 0) cout << "Complex Eigenproblem time: " << tm.real() << endl;

    }

//  Gram-Schmidt orthogonalization of matrix a
    for ( int m = 0; m < a.nblocks(); m++ )
      for ( int l = 0; l < a.mblocks(); l++ )
        for ( int y = 0; y < a.nbs(m); y++ )
          for ( int x = 0; x < a.mbs(l); x++ )
          {
            int i = a.i(l,x);
            int j = a.j(m,y);
            //double aij = aa(i,j);
            int iii = x + l*a.mb();
            int jjj = y + m*a.nb();
            int ival = iii + jjj * a.mloc();
            if ( i == j )
              a[ival] = i + 1.e-6*drand48();
            else
              a[ival] = 1.e-6*drand48();
          }
    DoubleMatrix s(a.context(),a.n(),a.n(),a.nb(),a.nb());
    tm.reset();
    tm.start();
    s.syrk('l','t',2.0,a,0.0);
    if (mype == 0) cout << "Gram syrk time: " << tm.real() << endl;

    tm.reset();
    tm.start();
    s.syr('l',-1.0,a,0,'r');
    tm.stop();
    if (mype == 0) cout << "Gram syr time: " << tm.real() << endl;

    // Cholesky decomposition
    tm.reset();
    tm.start();
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
    tm.stop();
    if (mype == 0) cout << "Gram Cholesky time: " << tm.real() << endl;

    // Triangular solve
    tm.reset();
    tm.start();
    // solve triangular system X * L^T = C
    a.trsm('r','l','t','n',1.0,s);
    tm.stop();
    if (mype == 0) cout << "Gram triangular solve time: " << tm.real() << endl;
  }

  MPI_Finalize();
}
