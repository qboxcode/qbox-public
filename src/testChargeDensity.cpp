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
// testChargeDensity.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "Context.h"
#include "Wavefunction.h"
#include "ChargeDensity.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Timer.h"

#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

#include <mpi.h>

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  {
    // use:
    // testChargeDensity a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel nkp nkspin
    //  ngb nstb nkpb nspb
    if ( argc != 18 )
    {
      cout << "use: testChargeDensity a0 a1 a2 b0 b1 b2 c0 c1 c2 "
           << "ecut nel nkp nspin ngb nstb nkpb nspb"
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
    int nkp = atoi(argv[12]);
    int nspin = atoi(argv[13]);
    int ngb = atoi(argv[14]);
    int nstb = atoi(argv[15]);
    int nkpb = atoi(argv[16]);
    int nspb = atoi(argv[17]);

    MPIdata::set(ngb,nstb,nkpb,nspb);
    cout << MPIdata::rank() << ": ngb=" << ngb << " nstb=" << nstb
         << " nkpb=" << nkpb << " nspb=" << nspb << endl;
    cout << MPIdata::rank() << ": igb=" << MPIdata::igb()
         << " istb=" << MPIdata::istb()
         << " ikpb=" << MPIdata::ikpb()
         << " ispb=" << MPIdata::ispb() << endl;

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

    for ( int ikp = 1; ikp < nkp; ikp++ )
    {
      wf.add_kpoint(D3vector((0.5*(ikp))/(nkp-1),0,0),1.0);
    }

    wf.info(cout,"wavefunction");

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

    wf.update_occ(0.0);

    // compute charge density in real space
    Timer tmrho;
    tmrho.reset();
    tmrho.start();
    cout << " ChargeDensity::ctor..." << endl;
    ChargeDensity cd(wf);
    tmrho.stop();
    cout << " ChargeDensity::ctor: CPU/Real: "
         << tmrho.cpu() << " / " << tmrho.real() << endl;

    tmrho.reset();
    tmrho.start();
    cout << " ChargeDensity::update_density..." << endl;
    cd.update_density();
    tmrho.stop();
    cout << " ChargeDensity::update_density: CPU/Real: "
         << tmrho.cpu() << " / " << tmrho.real() << endl;

    cout << cd;

    tmrho.reset();
    tmrho.start();
    cout << " ChargeDensity::update_rhor..." << endl;
    cd.update_rhor();
    tmrho.stop();
    cout << " ChargeDensity::update_rhor: CPU/Real: "
         << tmrho.cpu() << " / " << tmrho.real() << endl;

    cout << cd;
  }
  MPI_Finalize();
}
