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
// testWavefunction.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "Context.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Timer.h"

#include <iostream>
#include <cassert>
#include <cstdlib> // atoi
using namespace std;

#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  {
    // use:
    // testWavefunction a0 a1 a2 b0 b1 b2 c0 c1 c2
    //                  ecut nel nempty nspin nkp ngb nstb nspb nkpb
    if ( argc != 19 )
    {
      cout << "use: testWavefunction a0 a1 a2 b0 b1 b2 c0 c1 c2 "
           << "ecut nel nempty nspin nkp ngb nstb nspb nkpb"
      << endl;
      return 1;
    }
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    cout << " volume: " << cell.volume() << endl;
    double ecut = atof(argv[10]);
    int nel = atoi(argv[11]);
    int nempty = atoi(argv[12]);
    int nspin = atoi(argv[13]);
    int nkp = atoi(argv[14]);
    int ngb = atoi(argv[15]);
    int nstb = atoi(argv[16]);
    int nspb = atoi(argv[17]);
    int nkpb = atoi(argv[18]);

    MPIdata::set(ngb,nstb,nspb,nkpb);
    cout << MPIdata::rank() << ": ngb=" << ngb << " nstb=" << nstb
         << " nspb=" << nspb << " nkpb=" << nkpb << endl;
    cout << MPIdata::rank() << ": igb=" << MPIdata::igb()
         << " istb=" << MPIdata::istb()
         << " ispb=" << MPIdata::ispb()
         << " ikpb=" << MPIdata::ikpb() << endl;

    Context sd_ctxt(MPIdata::sd_comm(),ngb,nstb);
    Wavefunction wf(sd_ctxt);

    Timer tm;

    tm.reset(); tm.start();
    wf.resize(cell,cell,ecut);
    tm.stop();
    cout << " wf.resize: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset(); tm.start();
    wf.set_nel(nel);
    tm.stop();
    cout << " wf.set_nel: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset(); tm.start();
    wf.set_nspin(nspin);
    tm.stop();
    cout << " wf.set_nspin: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset(); tm.start();
    wf.set_nempty(nempty);
    tm.stop();
    cout << " wf.set_nempty: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    for ( int ikp = 0; ikp < nkp-1; ikp++ )
    {
      wf.add_kpoint(D3vector((0.5*(ikp+1))/(nkp-1),0,0),1.0);
    }

    tm.reset();
    tm.start();
    wf.randomize(0.1);
    tm.stop();
    cout << " wf.randomize: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset();
    tm.start();
    wf.gram();
    cout << " wf.gram: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    cout << " copy constructor...";
    Wavefunction wfm(wf);
    cout << "done" << endl;

    cout << " wfm gram ...";
    wfm.gram();
    cout << "done" << endl;

    cout << " randomize ...";
    wf.randomize(0.1);
    cout << "done" << endl;

    cout << " update_occ ...";
    wf.update_occ(0.0);
    cout << "done" << endl;

    wf.info(cout,"wavefunction");
  }
  MPI_Finalize();
}
