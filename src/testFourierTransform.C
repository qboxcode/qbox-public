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
// testFourierTransform.C
//
// test and timing of the Qbox FourierTransform class
//
// The FourierTransform functionality is tested in the following 8 tests
// relevant to the use of the class in different parts of Qbox

#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;

#include "Basis.h"
#include "FourierTransform.h"
#include "Timer.h"

#if USE_APC
#include "apc.h"
#endif

double fft_flops(int n)
{
  return 5.0 * n * log((double) n) / log(2.0);
}

int main(int argc, char **argv)
{
  Timer tm;
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
#if USE_APC
  ApcInit();
#endif
  {

#if USE_MPI
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);
#else
  const char* processor_name = "serial";
#endif

  int mype;
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

  // start scope of wf-v transforms
  {
  // transform and interpolate as for wavefunctions

  Basis basis(MPI_COMM_WORLD,kpoint);
  basis.resize(cell,cell,ecut);
  FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
  FourierTransform ft2(basis,2*basis.np(0),2*basis.np(1),2*basis.np(2));
  vector<complex<double> > f1(ft.np012loc());
  vector<complex<double> > f2(ft2.np012loc());
  vector<complex<double> > x(basis.localsize());
  vector<complex<double> > x1(basis.localsize());
  vector<complex<double> > x2(basis.localsize());

  double flops = 2*basis.nrod_loc() *      fft_flops(ft2.np2()) +
                 ft2.np1()/2 * ft2.np2() * fft_flops(ft2.np0()) +
                 ft2.np0()   * ft2.np2() * fft_flops(ft2.np1());
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
#if 1
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
#endif

  // test ft (small grid)
  cout << mype << ": ft.np2_loc(): " << ft.np2_loc() << endl;
  MPI_Barrier(MPI_COMM_WORLD);
  cout << " test ft: ";
  ft.forward(&f1[0],&x[0]);
  cout << " forward done ";
  ft.backward(&x[0],&f1[0]);
  cout << " backward done " << endl;
  MPI_Barrier(MPI_COMM_WORLD);

#if 1

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(1);
#endif
  ft2.forward(&f2[0],&x[0]);
#if USE_APC
  ApcStop(1);
#endif
  tm.stop();
  cout << " fwd1: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " fwd1: vgrid->wf" << endl;
  cout << " fwd1: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
  cout << " fwd1: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
  cout << " fwd1: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
  cout << " fwd1: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
  cout << " fwd1: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
  cout << " fwd1: tm_f_map:    " << ft2.tm_f_map.real() << endl;
  cout << " fwd1: tm_f_total:  " << ft2.tm_f_fft.real() +
                                    ft2.tm_f_mpi.real() +
                                    ft2.tm_f_pack.real() +
                                    ft2.tm_f_unpack.real() +
                                    ft2.tm_f_zero.real() +
                                    ft2.tm_f_map.real() << endl;
  cout << " fwd1 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(2);
#endif
  ft2.backward(&x[0],&f2[0]);
#if USE_APC
  ApcStop(2);
#endif
  tm.stop();
  cout << " bwd1: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " bwd1: wf->vgrid" << endl;
  cout << " bwd1: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
  cout << " bwd1: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
  cout << " bwd1: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
  cout << " bwd1: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
  cout << " bwd1: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
  cout << " bwd1: tm_b_map:    " << ft2.tm_b_map.real() << endl;
  cout << " bwd1: tm_b_total:  " << ft2.tm_b_fft.real() +
                                    ft2.tm_b_mpi.real() +
                                    ft2.tm_b_pack.real() +
                                    ft2.tm_b_unpack.real() +
                                    ft2.tm_b_zero.real() +
                                    ft2.tm_b_map.real() << endl;
  cout << " bwd1 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(3);
#endif
  ft2.forward(&f2[0],&x[0]);
#if USE_APC
  ApcStop(3);
#endif
  tm.stop();
  cout << " fwd2: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " fwd2: vgrid->wf" << endl;
  cout << " fwd2: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
  cout << " fwd2: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
  cout << " fwd2: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
  cout << " fwd2: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
  cout << " fwd2: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
  cout << " fwd2: tm_f_map:    " << ft2.tm_f_map.real() << endl;
  cout << " fwd2: tm_f_total:  " << ft2.tm_f_fft.real() +
                                    ft2.tm_f_mpi.real() +
                                    ft2.tm_f_pack.real() +
                                    ft2.tm_f_unpack.real() +
                                    ft2.tm_f_zero.real() +
                                    ft2.tm_f_map.real() << endl;

  //cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
  //     << " " << 2*basis.np(2) << " ";
  cout << " fwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(4);
#endif
  ft2.backward(&x[0],&f2[0]);
#if USE_APC
  ApcStop(4);
#endif
  tm.stop();
  cout << " bwd2: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " bwd2: wf->vgrid" << endl;
  cout << " bwd2: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
  cout << " bwd2: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
  cout << " bwd2: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
  cout << " bwd2: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
  cout << " bwd2: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
  cout << " bwd2: tm_b_map:    " << ft2.tm_b_map.real() << endl;
  cout << " bwd2: tm_b_total:  " << ft2.tm_b_fft.real() +
                                    ft2.tm_b_mpi.real() +
                                    ft2.tm_b_pack.real() +
                                    ft2.tm_b_unpack.real() +
                                    ft2.tm_b_zero.real() +
                                    ft2.tm_b_map.real() << endl;

  //cout << " " << 2*basis.np(0) << " " << 2*basis.np(1)
  //     << " " << 2*basis.np(2) << " ";
  cout << " bwd2 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  // double transform
  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(5);
#endif
  ft2.forward(&f2[0],&x1[0],&x2[0]);
#if USE_APC
  ApcStop(5);
#endif
  tm.stop();
  cout << " fwd3: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " fwd3: vgrid->wf double transform" << endl;
  cout << " fwd3: tm_f_fft:    " << ft2.tm_f_fft.real() << endl;
  cout << " fwd3: tm_f_mpi:    " << ft2.tm_f_mpi.real() << endl;
  cout << " fwd3: tm_f_pack:   " << ft2.tm_f_pack.real() << endl;
  cout << " fwd3: tm_f_unpack: " << ft2.tm_f_unpack.real() << endl;
  cout << " fwd3: tm_f_zero:   " << ft2.tm_f_zero.real() << endl;
  cout << " fwd3: tm_f_map:    " << ft2.tm_f_map.real() << endl;
  cout << " fwd3: tm_f_total:  " << ft2.tm_f_fft.real() +
                                    ft2.tm_f_mpi.real() +
                                    ft2.tm_f_pack.real() +
                                    ft2.tm_f_unpack.real() +
                                    ft2.tm_f_zero.real() +
                                    ft2.tm_f_map.real() << endl;
  cout << " fwd3 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  ft2.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(6);
#endif
  ft2.backward(&x1[0],&x2[0],&f2[0]);
#if USE_APC
  ApcStop(6);
#endif
  tm.stop();
  cout << " bwd3: size: " << ft2.np0() << " "
       << ft2.np1() << " " << ft2.np2() << endl;
  cout << " bwd3: wf->vgrid double transform" << endl;
  cout << " bwd3: tm_b_fft:    " << ft2.tm_b_fft.real() << endl;
  cout << " bwd3: tm_b_mpi:    " << ft2.tm_b_mpi.real() << endl;
  cout << " bwd3: tm_b_pack:   " << ft2.tm_b_pack.real() << endl;
  cout << " bwd3: tm_b_unpack: " << ft2.tm_b_unpack.real() << endl;
  cout << " bwd3: tm_b_zero:   " << ft2.tm_b_zero.real() << endl;
  cout << " bwd3: tm_b_map:    " << ft2.tm_b_map.real() << endl;
  cout << " bwd3: tm_b_total:  " << ft2.tm_b_fft.real() +
                                    ft2.tm_b_mpi.real() +
                                    ft2.tm_b_pack.real() +
                                    ft2.tm_b_unpack.real() +
                                    ft2.tm_b_zero.real() +
                                    ft2.tm_b_map.real() << endl;
  cout << " bwd3 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*flops/tm.real() << " MFlops" << endl;

#endif
  } // end of scope for wf-v transforms

#if 1
  ////////////////////////////////////////////////////////////
  // v(g)->vgrid
  Basis vbasis(MPI_COMM_WORLD,kpoint);
  vbasis.resize(cell,cell,4.0*ecut);
  cout << " vbasis.np() = " << vbasis.np(0) << " " << vbasis.np(1)
       << " " << vbasis.np(2) << endl;
  FourierTransform vft(vbasis,vbasis.np(0),vbasis.np(1),vbasis.np(2));
  vector<complex<double> > vf(vft.np012loc());
  vector<complex<double> > vg(vbasis.localsize());

  double vflops = 2*vbasis.nrod_loc() *      fft_flops(vft.np2()) +
                   vft.np1()   * vft.np2() * fft_flops(vft.np0()) +
                   vft.np0()   * vft.np2() * fft_flops(vft.np1());
  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  vft.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(7);
#endif
  vft.forward(&vf[0],&vg[0]);
#if USE_APC
  ApcStop(7);
#endif
  tm.stop();
  cout << " fwd4: size: " << vft.np0() << " "
       << vft.np1() << " " << vft.np2() << endl;
  cout << " fwd4: vgrid->v(g)" << endl;
  cout << " fwd4: tm_f_fft:    " << vft.tm_f_fft.real() << endl;
  cout << " fwd4: tm_f_mpi:    " << vft.tm_f_mpi.real() << endl;
  cout << " fwd4: tm_f_pack:   " << vft.tm_f_pack.real() << endl;
  cout << " fwd4: tm_f_unpack: " << vft.tm_f_unpack.real() << endl;
  cout << " fwd4: tm_f_zero:   " << vft.tm_f_zero.real() << endl;
  cout << " fwd4: tm_f_map:    " << vft.tm_f_map.real() << endl;
  cout << " fwd4: tm_f_total:  " << vft.tm_f_fft.real() +
                                    vft.tm_f_mpi.real() +
                                    vft.tm_f_pack.real() +
                                    vft.tm_f_unpack.real() +
                                    vft.tm_f_zero.real() +
                                    vft.tm_f_map.real() << endl;
  cout << " fwd4 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*vflops/tm.real() << " MFlops" << endl;

  MPI_Barrier(MPI_COMM_WORLD);
  tm.reset();
  vft.reset_timers();
  tm.start();
#if USE_APC
  ApcStart(8);
#endif
  vft.backward(&vg[0],&vf[0]);
#if USE_APC
  ApcStop(8);
#endif
  tm.stop();
  cout << " bwd4: size: " << vft.np0() << " "
       << vft.np1() << " " << vft.np2() << endl;
  cout << " bwd4: v(g)->vgrid" << endl;
  cout << " bwd4: tm_b_fft:    " << vft.tm_b_fft.real() << endl;
  cout << " bwd4: tm_b_mpi:    " << vft.tm_b_mpi.real() << endl;
  cout << " bwd4: tm_b_pack:   " << vft.tm_b_pack.real() << endl;
  cout << " bwd4: tm_b_unpack: " << vft.tm_b_unpack.real() << endl;
  cout << " bwd4: tm_b_zero:   " << vft.tm_b_zero.real() << endl;
  cout << " bwd4: tm_b_map:    " << vft.tm_b_map.real() << endl;
  cout << " bwd4: tm_b_total:  " << vft.tm_b_fft.real() +
                                    vft.tm_b_mpi.real() +
                                    vft.tm_b_pack.real() +
                                    vft.tm_b_unpack.real() +
                                    vft.tm_b_zero.real() +
                                    vft.tm_b_map.real() << endl;
  cout << " bwd4 time: " << tm.cpu() << " / " << tm.real()
  << "    " << 1.e-6*vflops/tm.real() << " MFlops" << endl;
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
#endif

  } // Context scope
#if USE_APC
  ApcFinalize();
#endif
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
