////////////////////////////////////////////////////////////////////////////////
//
// testContext.c
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testContext.C,v 1.4 2003-01-10 00:29:38 fgygi Exp $

#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#ifdef USE_MPI  
#include <mpi.h>
#endif

#include "Context.h"

int main(int argc, char **argv)
{
  int mype;
  int npes;
#ifdef USE_MPI  

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  { // start Context scope
  
  // test reference counting
  if ( npes%2 == 0 && npes >= 4 )
  {
    Context *c1 = new Context(npes/2,2);
    cout << "c1: " << *c1 << endl;
    Context *c2;
    if ( c1->active() )
      c2 = new Context(*c1,npes/2,npes/2,1);
    else
      c2 = 0;
    // this line causes crash: Context *c2 = new Context(*c1,1,1,1);
    delete c1;
    if ( c2 != 0 ) cout << c2->mype() << " c2: " << *c2;
    delete c2;
  }
  
#if 0
  Context c_global;
  
  if ( npes%2 == 0 && npes >= 2 )
  {
    Context sub(c_global,c_global.size()/2,2);
    Context half0(sub,0,c_global.size()/2,1);
    Context half1(sub,c_global.size()/2,c_global.size()/2,1);
    cout << c_global.mype() << ": half0: " << half0;
    cout << c_global.mype() << ": half1: " << half1;
  }
  
  if ( npes >= 16 )
  {
    Context subcontext(c_global,npes/4,4,'c');
    cout << c_global.mype() << ": subcontext: " << subcontext;
    
    if ( subcontext.active() )
    {
      Context sd0(subcontext,0,0,4,2,4,2);
      Context sd1(subcontext,0,2,4,2,4,2);
      cout << c_global.mype() << ": sd0: " << sd0;
      cout << c_global.mype() << ": sd1: " << sd1;
    }
  }
    
  if ( npes == 6 )
  {
  
    Context rowmajor(2,3,'R');
    cout << " rowmajor: " << rowmajor;
    Context colmajor(2,3,'C');
    cout << " colmajor: " << colmajor;
    cout << " colmajor.mype(): " << colmajor.mype()
         << " myproc(): " << colmajor.myproc() 
         << " myrow(): " << colmajor.myrow() 
         << " mycol(): " << colmajor.mycol() 
         << endl;
    
    vector<Context*> c;
    
    c.push_back( new Context(c_global,3,2,'r') );
    c.push_back( new Context(c_global,'s') );
    c.push_back( new Context(c_global,1,3,1,'c') );
    c.push_back( new Context(*c[0],1,0,2,2,4,1) );
    c.push_back( new Context(*c[0],'r',c[0]->myrow()) );
    c.push_back( new Context(*c[0],'c',c[0]->mycol()) );
    c.push_back( new Context() );
    
    for ( int i = 0; i < c.size(); i++ )
    {
      if ( c[i]->myproc() == 0 )
        cout << "c[" << i << "]: " << *c[i];
    }
    
    Context c1(2,3,'c');
    Context c2 = c1;
    Context c3(c1);
    Context c4;
    c4 = c1;
    
    cout << " c1: " << c1;
    cout << " c2: " << c2;
    cout << " c3: " << c3;
    cout << " c4: " << c4;
    
    vector<Context> cv;
    cv.resize(4);
    cv.push_back(c1);
    cv.push_back(c2);
    cv.push_back(c3);
    cv.push_back(c4);
    
    for ( int i = 0; i < cv.size(); i++ )
      cout << " cv[" << i << "]: " << cv[i];

    // test dgsum2d function
    double a = c[0]->mype();
    cout << c[0]->mype() << ": a     = " << a << endl;
    c[0]->dgsum2d('R',1,1,&a,1);
    cout << c[0]->mype() << ": a_sum_row = " << a << endl;
    c[0]->dgsum2d('C',1,1,&a,1);
    cout << c[0]->mype() << ": a_sum_all = " << a << endl;
    
    
  }
#endif
  
  } // end Context scope

  MPI_Finalize();

#else

  mype=0;
  npes=1;
  {
    BlacsContext c1(1,1);
    cout << " c1.ictxt = " << c1.ictxt() << endl;
  }

#endif
}
