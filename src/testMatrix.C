// $Id: testMatrix.C,v 1.8 2004-11-30 23:00:09 fgygi Exp $
//
// test Matrix
//
// multiply a matrix a(m,k) by b(k,n) to get c(m,n)
// using blocks of size (mb,nb) on a process grid (nprow,npcol)
//
// use: testMatrix nprow npcol input_file [-check] 
// input_file:
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

#include "Timer.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "Context.h"
#include "Matrix.h"

int aa(int i, int j) { return i+2*j; }
int bb(int i, int j) { return i-j-3; }

int main(int argc, char **argv)
{
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

  char* infilename = argv[3];
  ifstream infile(infilename);

  assert(argc == 4 || argc == 5);
  bool tcheck = false;
  if ( argc == 5 )
  {
    if ( !strcmp(argv[4],"-check") )
      tcheck = true;
    else
    {
      cerr << " invalid argv[4]" << endl;
#if USE_MPI
      MPI_Abort(MPI_COMM_WORLD,2);
#else
      exit(2);
#endif
    }
  }
      
  int nprow = atoi(argv[1]);
  int npcol = atoi(argv[2]);
  
  Timer tm;
  
  int m_a, n_a, mb_a, nb_a;
  int m_b, n_b, mb_b, nb_b;
  int m_c, n_c, mb_c, nb_c;
  char ta, tb;
  if(mype == 0)
  {
    infile >> m_a >> n_a >> mb_a >> nb_a >> ta;
    cout<<"m_a="<<m_a<<", n_a="<<n_a<<endl;
    infile >> m_b >> n_b >> mb_b >> nb_b >> tb;
    cout<<"m_b="<<m_b<<", n_b="<<n_a<<endl;
    infile >> m_c >> n_c >> mb_c >> nb_c;
    cout<<"m_c="<<m_c<<", n_c="<<n_c<<endl;
  }
#ifdef USE_MPI
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
#endif
  {  
    if ( ta == 'N' ) ta = 'n';
    if ( tb == 'N' ) tb = 'n';

    Context ctxt(nprow,npcol);

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
    c.gemm(ta,tb,1.0,a,b,0.0);
    tm.stop();


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

    cout << " CPU/Real: " << setw(8) << tm.cpu() 
         << " / " << setw(8) << tm.real();
    if ( tm.real() > 0.0 )
    {
      int kmax = ( ta == 'n' ) ? a.n() : a.m();
      cout << "  MFlops: " 
           << (2.0e-6*m_c*n_c*kmax) / tm.real() << endl;
    }
    
    double norma=a.nrm2();
    if(mype == 0)cout<<"Norm(a)="<<norma<<endl;
    if(mype == 0)cout<<"DoubleMatrix::matgather..."<<endl;
    double*  aa=new double[a.m()*a.n()];
    a.matgather(aa, a.m());
    if(mype == 0)cout<<"DoubleMatrix::init..."<<endl;
    b.init(aa, a.m());
    double norm=b.nrm2();
    if ( mype == 0 ) cout << "Norm(b)=" << norm << endl;
    if ( fabs(norm-norma)>0.000001 )
       cout << "DoubleMatrix: problem with matgather/init" << endl;

    if ( c.n() == b.m() && c.m() == b.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::transpose..."<<endl;
      c.transpose(1.0,b,0.0);
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
    
    if ( a.m() == c.m() && a.n() == c.n() )
    {
      if(mype == 0)cout<<"DoubleMatrix::operator=..."<<endl;
      c=a;
    }
    
    if(mype == 0)cout<<"DoubleMatrix::nrm2..."<<endl;
    norm=c.nrm2();
    if (mype == 0) cout<<"Norm="<<norm<<endl;
    
    if ( mype == 0 ) cout << "DoubleMatrix::syrk...";
    ComplexMatrix d(ctxt,m_c,n_c,mb_c,nb_c);
    DoubleMatrix e(d);
    c.syrk('l','t',2.0,e,0.0);
    if ( mype == 0 ) cout << "done" << endl;
    
    a.identity();
    DoubleMatrix a2(a);
    a -= a2;
    norm = a.nrm2();
    if (mype == 0) cout << "Norm(a)=" << norm << endl;

    if ( ctxt.nprow() >= 2 )
    {
      Context subctxt(2,1);
      DoubleMatrix csub(subctxt);
      csub.resize(c.m(),c.n());
      csub.getsub(c,c.m(),c.n(),0,0);
      norm = c.nrm2();
      if (mype == 0) cout << "Norm(c)=" << norm << endl;
      norm = csub.nrm2();
      if (mype == 0) cout << "Norm(csub)=" << norm << endl;
    }
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
}
