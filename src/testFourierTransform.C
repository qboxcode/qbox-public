//
// testFourierTransform.C
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "Context.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "Timer.h"

int main(int argc, char **argv)
{
  Timer tm;
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that Context objects get destructed before
  // the MPI_Finalize call
  {
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);

  Context ctxt_global;
  int mype = ctxt_global.mype();
  //cout << " Process " << mype << " on " << processor_name << endl;

  if ( argc != 14 )
  {
    cout << 
    " use: testFourierTransform "
    << "a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut kx ky kz" << endl;
  }
  D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
  D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
  D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
  UnitCell cell(a,b,c);
  const double omega = cell.volume();
  double ecut = atof(argv[10]);
  D3vector kpoint(atof(argv[11]),atof(argv[12]),atof(argv[13]));
  
  //cout << " ctxt_global: " << ctxt_global;
  //cout << " ctxt_global.comm(): " << ctxt_global.comm() << endl;
  Context ctxt(ctxt_global,'c');
  //cout << " ctxt: " << ctxt;
  //cout << " ctxt.comm(): " << ctxt.comm() << endl;
  Basis basis(ctxt,kpoint);
  basis.resize(cell,cell,ecut);
  
  FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
  
  vector<complex<double> > x(basis.localsize());
  double flops = 5.0 * ft.np012() * 
    ( log((double)basis.np(0)) + log((double)basis.np(1)) + 
      log((double)basis.np(2)) ) / log(2.0);
  if ( ctxt.onpe0() )
  {
    cout << " basis.localsize() = " << basis.localsize() << endl;
    cout << " basis.np() = " << basis.np(0) << " " << basis.np(1)
         << " " << basis.np(2) << endl;
    cout << " flop count: " << flops << endl;
  }
  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::right,ios::adjustfield);
  cout << setprecision(6);
  
  const double rc = 1.0;
#if 1
  // Initialize with Fourier coefficients of a normalized gaussian distribution
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    double y = 1.0/omega * exp( -0.25 * g2 * rc*rc );
    x[i] = y;
    // x[i] = complex<double>(y,y);
  }
#endif

#if 0
  x[0] = 1.0;
#endif
  
  vector<complex<double> > f(ft.np012loc());
  
  tm.reset();
  tm.start();
  ft.backward(&x[0],&f[0]);
  tm.stop();
  cout << " " << basis.np(0) << " " << basis.np(1)
       << " " << basis.np(2) << " ";
  cout << " bwd time: " << tm.cpu() << " / " << tm.real() 
       << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  
#if 0
  for ( int i = 0; i < basis.localsize(); i++ )
    cout << setw(8) << basis.gx(i) << " " 
         << setw(8) << basis.gx(i+basis.localsize()) << " " 
         << setw(8) << basis.gx(i+2*basis.localsize())
         << "     " << x[i] << endl;
  for ( int i = 0; i < ft.np0(); i++ )
    for ( int j = 0; j < ft.np1(); j++ )
      for ( int k = 0; k < ft.np2_loc(); k++ )
        cout << mype << ": "
             << i << " " << j << " " << k+ft.np2_first() << " "
             << f[ft.index(i,j,k)] << endl;
#endif
             
  // integral in r space must be 1.0
  double sum,tsum = 0.0;
  for ( int i = 0; i < f.size(); i++ ) 
    tsum += real(f[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  if ( ctxt.onpe0() )
    cout << " sum: " << sum * omega / ft.np012() << endl;
  
  FourierTransform ft2(basis,2*basis.np(0),2*basis.np(1),2*basis.np(2));
  vector<complex<double> > f2(ft2.np012loc());
  tm.reset();
  tm.start();
  ft2.backward(&x[0],&f2[0]);
  tm.stop();
  
  // integral in r space must be 1.0
  tsum = 0.0;
  for ( int i = 0; i < f2.size(); i++ ) 
    tsum += real(f2[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  

  cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
       << " " << 2*basis.np(2) << " ";
  cout << " bwd time: " << tm.cpu() << " / " << tm.real() 
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  cout << " sum on f2: " << sum * omega / ft2.np012() << endl;
  
#if 0
  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized plane wave
  //////////////////////////////////////////////////////////////////////////////
  
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    x[i] = 0.0;
  }
  if ( ctxt.myproc() == 0 ) x[1] = 1.0/sqrt(2.0);
  
  ft.backward(&x[0],&f[0]);

#if 0
  for ( int i = 0; i < basis.localsize(); i++ )
    cout << basis.kv(3*i) << " " << basis.kv(3*i+1) << " " << basis.kv(3*i+2)
         << "     " << x[i] << endl;
  for ( int i = 0; i < ft.np0(); i++ )
    for ( int j = 0; j < ft.np1(); j++ )
      for ( int k = 0; k < ft.np2_loc(); k++ )
        cout << mype << ": "
             << i << " " << j << " " << k+ft.np2_first() << " "
             << f[ft.index(i,j,k)] << endl;
#endif
             
  // integral of f^2 in r space must be 1.0
  tsum = 0.0;
  for ( int i = 0; i < f.size(); i++ ) 
    tsum += norm(f[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  cout << " sum pw^2: " << sum / ft.np012() << endl;
  
  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized gaussian
  //////////////////////////////////////////////////////////////////////////////
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    x[i] = 1.0 / sqrt(omega) * pow(2.0*M_PI*rc*rc,0.75) * 
           exp( -0.25 * g2 * rc*rc );
  }
  
  ft2.backward(&x[0],&f2[0]);
  
//   for ( int i = 0; i < basis.localsize(); i++ )
//     cout << basis.kv(3*i) << " " << basis.kv(3*i+1) << " " << basis.kv(3*i+2)
//          << "     " << x[i] << endl;
//   for ( int i = 0; i < ft2.np0(); i++ )
//     for ( int j = 0; j < ft2.np1(); j++ )
//       for ( int k = 0; k < ft2.np2_loc(); k++ )
//         cout << mype << ": "
//              << i << " " << j << " " << k+ft2.np2_first() << " "
//              << f2[ft2.index(i,j,k)] << endl;
             
  // integral of gaussian^2 in r space must be 1.0
  tsum = 0.0;
  for ( int i = 0; i < f.size(); i++ ) 
    tsum += norm(f[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());
  
  cout << " sum gaussian^2: " << sum / ft.np012() << endl;
  
#endif

  //////////////////////////////////////////////////////////////////////////////
  // Test forward transform
  //////////////////////////////////////////////////////////////////////////////
  for ( int i = 0; i < ft2.np0(); i++ )
    for ( int j = 0; j < ft2.np1(); j++ )
      for ( int k = 0; k < ft2.np2_loc(); k++ )
         f2[ft2.index(i,j,k)] = cos(k*2*M_PI/ft2.np2()) + sin(i*2*M_PI/ft2.np0()) ;
  ft2.forward(&f2[0],&x[0]);
#if 0
  for ( int i = 0; i < basis.localsize(); i++ )
    cout << basis.kv(3*i) << " " << basis.kv(3*i+1) << " " << basis.kv(3*i+2)
         << "   x[G]= " << x[i] << endl;
#endif
  tm.reset();
  tm.start();
  ft2.forward(&f2[0],&x[0]);
  tm.stop();
  cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
       << " " << 2*basis.np(2) << " ";
  cout << " fwd time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  tm.reset();
  tm.start();
  ft2.backward(&x[0],&f2[0]);
  tm.stop();
  cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
       << " " << 2*basis.np(2) << " ";
  cout << " bwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
  
  tm.reset();
  tm.start();
  ft2.forward(&f2[0],&x[0]);
  tm.stop();
  cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
       << " " << 2*basis.np(2) << " ";
  cout << " fwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
