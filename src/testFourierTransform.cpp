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
//
// testFourierTransform.cpp
//
// test and timing of the Qbox FourierTransform class
//
// The FourierTransform functionality is tested in the following tests
// relevant to the use of the class in different parts of Qbox

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
using namespace std;

#include "Basis.h"
#include "FourierTransform.h"
#include "Timer.h"

int mype, npes;

double fft_flops_1d(int n)
{
  return 5.0 * n * log((double) n) / log(2.0);
}

double fft_flops(Basis& b, FourierTransform& ft)
{
  return  2*b.nrod_loc() * fft_flops_1d(ft.np2()) +
          ft.np1()/2 * ft.np2() * fft_flops_1d(ft.np0()) +
          ft.np0()   * ft.np2() * fft_flops_1d(ft.np1());
}

void print_timing(std::string name, FourierTransform& ft,
  double flops, Timer& tm)
{
  if ( mype == 0 )
  {
    cout << "*****************************************************" << endl;
    cout << " " << name << endl;
    cout << " size: " << ft.np0() << " " << ft.np1() << " " << ft.np2() << endl;
  }
  for ( int ipe = 0; ipe < npes; ipe++ )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if ( ipe == mype )
    {
      cout << ipe << ": tm_fwd:       " << ft.tm_fwd.real() << endl;
      cout << ipe << ": tm_bwd:       " << ft.tm_bwd.real() << endl;
      cout << ipe << ": tm_map_fwd:   " << ft.tm_map_fwd.real() << endl;
      cout << ipe << ": tm_map_bwd:   " << ft.tm_map_bwd.real() << endl;
      cout << ipe << ": tm_trans_fwd: " << ft.tm_trans_fwd.real() << endl;
      cout << ipe << ": tm_trans_bwd: " << ft.tm_trans_bwd.real() << endl;
      cout << ipe << ": tm_fxy:       " << ft.tm_fxy.real() << endl;
      cout << ipe << ": tm_fxy_inv:   " << ft.tm_fxy_inv.real() << endl;
      cout << ipe << ": tm_fz:        " << ft.tm_fz.real() << endl;
      cout << ipe << ": tm_fz_inv:    " << ft.tm_fz_inv.real() << endl;
      cout << ipe << ": time: " << tm.cpu() << " / " << tm.real()
      << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;
    }
  }
}

int main(int argc, char **argv)
{
  Timer tm;
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {

#if USE_MPI
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);
#else
  const char* processor_name = "serial";
#endif

  MPI_Comm_size(MPI_COMM_WORLD,&npes);
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  cout << " Process " << mype << " on " << processor_name << endl;

  D3vector a,b,c,kpoint;
  double ecut;
  if ( argc == 5 )
  {
    a = D3vector(atof(argv[1]),0,0);
    b = D3vector(0,atof(argv[2]),0);
    c = D3vector(0,0,atof(argv[3]));
    ecut = atof(argv[4]);
  }
  else if ( argc == 11 )
  {
    a = D3vector(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    b = D3vector(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    c = D3vector(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    ecut = atof(argv[10]);
  }
  else if ( argc == 14 )
  {
    a = D3vector(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    b = D3vector(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    c = D3vector(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    ecut = atof(argv[10]);
    kpoint = D3vector(atof(argv[11]),atof(argv[12]),atof(argv[13]));
  }
  else
  {
    cout << " use: testFourierTransform a b c ecut(a.u.) [kpoint] " << endl;
    return 1;
  }
  UnitCell cell(a,b,c);
  const double omega = cell.volume();

  // start scope of transforms
  {

  Basis basis(MPI_COMM_WORLD,kpoint);
  basis.resize(cell,cell,ecut);
  FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
  FourierTransform ft2(basis,2*basis.np(0),2*basis.np(1),2*basis.np(2));
  vector<complex<double> > f1(ft.np012loc());
  vector<complex<double> > f2(ft2.np012loc());
  vector<complex<double> > x(basis.localsize());
  vector<complex<double> > x1(basis.localsize());
  vector<complex<double> > x2(basis.localsize());

  double flops;
  if ( mype == 0 )
  {
    cout << " wfbasis.size() = " << basis.size() << endl;
    cout << " wfbasis.np() = " << basis.np(0) << " " << basis.np(1)
         << " " << basis.np(2) << endl;
    cout << " flop count: " << flops << endl;
  }
  cout << " wfbasis.nrod_loc(): " << basis.nrod_loc() << endl;
  cout << " zvec.size: "
       << 2*basis.nrod_loc()*ft2.np2() * sizeof(complex<double>)
       << endl;

  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::right,ios::adjustfield);
  cout << setprecision(6);

  const double rc = 1.0;
  // Initialize with Fourier coefficients of a normalized gaussian distribution
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    double y = 1.0/omega * exp( -0.25 * g2 * rc*rc );
    x[i] = y;
    x1[i] = y;
    x2[i] = y;
    // x[i] = complex<double>(y,y);
  }

  // wfgrid->wf transforms
  flops = fft_flops(basis,ft);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft.reset_timers();
  tm.start();
  ft.forward(&f1[0],&x[0]);
  tm.stop();
  print_timing("wfgrid->wf",ft,flops,tm);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft.reset_timers();
  tm.start();
  ft.backward(&x[0],&f2[0]);
  tm.stop();
  print_timing("wf->wfgrid",ft,flops,tm);

  // vgrid->wf transforms
  flops = fft_flops(basis,ft2);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.forward(&f2[0],&x[0]);
  tm.stop();
  print_timing("vgrid->wf",ft2,flops,tm);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.backward(&x[0],&f2[0]);
  tm.stop();
  print_timing("wf->vgrid",ft2,flops,tm);

  // double vgrid-wf transforms
  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.forward(&f2[0],&x1[0],&x2[0]);
  tm.stop();
  print_timing("vgrid->wf double transform",ft2,flops,tm);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
  ft2.backward(&x1[0],&x2[0],&f2[0]);
  tm.stop();
  print_timing("wf->vgrid double transform",ft2,flops,tm);

  // v(g)->vgrid transforms
  Basis vbasis(MPI_COMM_WORLD,kpoint);
  vbasis.resize(cell,cell,4.0*ecut);
  cout << " vbasis.np() = " << vbasis.np(0) << " " << vbasis.np(1)
       << " " << vbasis.np(2) << endl;
  FourierTransform vft(vbasis,vbasis.np(0),vbasis.np(1),vbasis.np(2));
  vector<complex<double> > vf(vft.np012loc());
  vector<complex<double> > vg(vbasis.localsize());

  flops = fft_flops(vbasis,vft);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  vft.reset_timers();
  tm.start();
  vft.forward(&vf[0],&vg[0]);
  tm.stop();
  print_timing("vgrid->vg",vft,flops,tm);

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  vft.reset_timers();
  tm.start();
  vft.backward(&vg[0],&vf[0]);
  tm.stop();
  print_timing("vg->vgrid",vft,flops,tm);

  } // end of scope for wf-v transforms

#if 0
  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized plane wave
  //////////////////////////////////////////////////////////////////////////////

  for ( int i = 0; i < basis.localsize(); i++ )
  {
    x[i] = 0.0;
  }
  if ( ctxt.myproc() == 0 ) x[1] = 1.0/sqrt(2.0);

  ft2.backward(&x[0],&f2[0]);

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

#if 0
  // integral of f^2 in r space must be 1.0
  double sum=0.0, tsum = 0.0;
  for ( int i = 0; i < f2.size(); i++ )
    tsum += norm(f2[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());

  cout << " sum pw^2: " << sum / ft2.np012() << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Integration of a 2-norm normalized gaussian
  //////////////////////////////////////////////////////////////////////////////
  for ( int i = 0; i < basis.localsize(); i++ )
  {
    double g2 = basis.g2(i);
    x[i] = 1.0 / sqrt(omega) * pow(2.0*M_PI*rc*rc,0.75) *
           exp( -0.25 * g2 * rc*rc );
  }
#endif

#if 0
  // Compute norm in g space
  double gnorm = 0.0;
  for ( int i = 0; i < basis.localsize(); i++ )
    gnorm += 2.0 * norm(x[i]);
  if ( ctxt.onpe0() )
    gnorm -= norm(x[0]);
  ctxt.dsum(1,1,&gnorm,1);
  cout << " gaussian gnorm: " << gnorm << endl;

  ft2.backward(&x[0],&f2[0]);
#endif

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
  for ( int i = 0; i < f2.size(); i++ )
    tsum += norm(f2[i]);
  MPI_Allreduce(&tsum,&sum,1,MPI_DOUBLE,MPI_SUM,ctxt.comm());

  cout << " gaussian rnorm: " << sum / ft2.np012() << endl;
#endif

  } // Context scope
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
