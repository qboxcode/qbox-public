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
// testEnergyFunctional.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Context.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "ChargeDensity.h"
#include "EnergyFunctional.h"
#include "Timer.h"

#include <iostream>
#include <iomanip>
#include <cassert>
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  {
    // use: testEnergyFunctional a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel
    if ( argc != 12 )
    {
      cout << "use: testEnergyFunctional a0 a1 a2 b0 b1 b2 c0 c1 c2 ecut nel"
           << endl;
      return 1;
    }
    D3vector a(atof(argv[1]),atof(argv[2]),atof(argv[3]));
    D3vector b(atof(argv[4]),atof(argv[5]),atof(argv[6]));
    D3vector c(atof(argv[7]),atof(argv[8]),atof(argv[9]));
    UnitCell cell(a,b,c);
    double ecut = atof(argv[10]);
    int nel = atoi(argv[11]);

    Timer tm;

    Context ctxt(MPI_COMM_WORLD);
    cout << " initial context: " << ctxt;
    Sample s(ctxt);

    s.wf.resize(cell,cell,ecut);
    s.wf.set_nel(nel);

    if ( ctxt.onpe0() ) cout << " nel: " << s.wf.nel() << endl;

    s.wf.update_occ(0.0);
    s.wf.randomize(0.05);

    tm.reset();
    tm.start();
    s.wf.gram();
    tm.stop();
    cout << " Gram: CPU/Real: " << tm.cpu() << " / " << tm.real() << endl;

    ChargeDensity cd(s.wf);
    tm.reset();
    tm.start();
    cout << " ChargeDensity::update_density..." << endl;
    cd.update_density();
    tm.stop();
    cout << " ChargeDensity::update_density: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    s.ctrl.xc = "LDA";
    s.ctrl.polarization = "OFF";
    tm.reset();
    tm.start();
    EnergyFunctional ef(s,cd);
    tm.stop();
    cout << " EnergyFunctional:ctor: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;

    tm.reset();
    tm.start();
    Wavefunction dwf(s.wf);
    vector<vector<double> > fion;
    valarray<double> sigma(6);
    double e = ef.energy(true,dwf,false,fion,false,sigma);
    cout << " ef.energy(): " << e << endl;
    tm.stop();
    cout << " EnergyFunctional:energy: CPU/Real: "
         << tm.cpu() << " / " << tm.real() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
}
