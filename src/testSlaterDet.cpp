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
// testSlaterDet.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Context.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Timer.h"
#include "MPIdata.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cstdlib> // atoi
using namespace std;

#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  int npes;
  MPI_Comm_size(MPI_COMM_WORLD,&npes);

  // use: testSlaterDet a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nst kx ky kz npr npc
  if ( argc != 17 )
  {
    cout <<
    "use: testSlaterDet a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nst kx ky kz npr npc"
    << endl;
    return 1;
  }
  double err;

  D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
  D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
  D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
  UnitCell cell(a,b,c);
  double ecut = atof(argv[10]);
  int nst = atoi(argv[11]);
  D3vector kpoint(atof(argv[12]),atof(argv[13]),atof(argv[14]));

  int npr = atoi(argv[15]);
  int npc = atoi(argv[16]);
  MPIdata::set(npr,npc);

  if ( MPIdata::onpe0() )
  {
    cout << " testSlaterDet" << endl;
    cout << " args: " << a << " " << b << " " << c << " "
         << ecut << " " << nst << " "
         << kpoint << " " << npr << " " << npc << endl;
  }

  std::map<std::string,Timer> tmap;

  { // Context scope
  Context ctxt(MPIdata::sd_comm(),npr,npc);
  if ( ctxt.myproc() == 0 )
    cout << " npr=" << npr << " npc=" << npc << endl;

  SlaterDet sd(ctxt,kpoint);

  sd.resize(cell,cell,ecut,nst);
  if ( ctxt.myproc() == 0 )
  {
    cout << " kpoint: " << sd.basis().kpoint() << endl;
    cout << " ngw:    " << sd.basis().size() << endl;
    cout << " nst:    " << sd.nst() << endl;
  }

  err = sd.ortho_error();
  if ( ctxt.myproc() == 0 )
    cout << " Initial orthogonality error before rand " << err << endl;

  sd.randomize(0.2/sd.basis().size());
  if ( ctxt.myproc() == 0 )
    cout << " sd.randomize: done" << endl;
  err = sd.ortho_error();
  if ( ctxt.myproc() == 0 )
    cout << " Orthogonality error after rand " << err << endl;

  tmap["gram"].reset();
  tmap["gram"].start();
  sd.gram();
  tmap["gram"].stop();

  err = sd.ortho_error();
  if ( ctxt.myproc() == 0 )
    cout << " Gram orthogonality error " << err << endl;

  SlaterDet sdm(ctxt,kpoint);
  sdm.resize(cell,cell,ecut,nst);

  sdm.c() = sd.c();
  sd.randomize(0.1/sd.basis().size());

  tmap["riccati"].reset();
  tmap["riccati"].start();
  sd.riccati(sdm);
  tmap["riccati"].stop();
  err = sd.ortho_error();
  if ( ctxt.myproc() == 0 )
    cout << " Riccati orthogonality error " << err << endl;

  sd.riccati(sdm);
  err = sd.ortho_error();
  if ( ctxt.myproc() == 0 )
  {
    cout << " Riccati orthogonality error " << err << endl;
    cout << " sd.total size (MB): " << setprecision(4)
         << sd.memsize() / 1048576.0 << endl;
    cout << " sd.local size (MB): " << setprecision(4)
         << sd.localmemsize() / 1048576.0 << endl;
  }

  // compute charge density in real space
  FourierTransform ft(sd.basis(),
    2*sd.basis().np(0), 2*sd.basis().np(1), 2*sd.basis().np(2));

  vector<complex<double> > f(ft.np012loc());
  vector<double> rho(ft.np012loc());

  if ( ctxt.myproc() == 0 )
    cout << " density" << endl;

  sd.update_occ(2*nst,1);
  double weight = 1.0;

  tmap["compute_density"].reset();
  tmap["compute_density"].start();
  sd.compute_density(ft,weight,&rho[0]);
  tmap["compute_density"].stop();

  // integral of rho in r space
  double sum = 0.0;
  for ( int i = 0; i < rho.size(); i++ )
    sum += rho[i];

  double tsum;
  MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,sd.context().comm());
  sum = tsum;
  if ( ctxt.myproc() == 0 )
    cout << " density sum: " << sum * cell.volume() / ft.np012() << endl;

  SlaterDet sdp(sd);
  //SlaterDet sdp(ctxt,kpoint);
  //sdp.resize(cell,cell,ecut,nst);
  //sdp = sd;
  sdp.randomize(0.001);
  sdp.gram();
  err = sdp.ortho_error();
  if ( ctxt.myproc() == 0 )
    cout << " Gram orthogonality error " << err << endl;

  Timer tmv;
  tmap["rs_mul_add"].reset();
  tmap["rs_mul_add"].start();
  sd.rs_mul_add(ft,&rho[0],sdp);
  tmap["rs_mul_add"].stop();

  //cout << sd;
#if 0
  ofstream outfile("sd.dat");
  tm.reset();
  tm.start();
  outfile << sd;
  tm.stop();
  cout << " write: CPU/Real: "
       << tm.cpu() << " / " << tm.real() << endl;
#endif

  ctxt.barrier();
  } // Context scope

  if ( MPIdata::onpe0() )
  {
    printf("\n");
    printf(" %-22s %12s %12s\n","timer name","min","max");
    printf(" %-22s %12s %12s\n","----------","---","---");
  }
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = i->second.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( MPIdata::onpe0() )
    {
      string s = (*i).first;
      printf(" %-22s %12.6f %12.6f\n",s.c_str(), tmin, tmax);
    }
  }
  MPI_Finalize();
}
