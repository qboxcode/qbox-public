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
// testSample.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: testSample.C,v 1.6 2009-11-30 02:23:26 fgygi Exp $

#include <iostream>
using namespace std;

#include "Context.h"
#include "SlaterDet.h"
#include "UnitCell.h"
#include "Sample.h"
#include "D3vector.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  {
    Context ctxt;

#if USE_MPI
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << ctxt.mype() << " on " << processor_name << endl;
#endif

    Sample s(ctxt);

    D3vector a(18, 0, 0);
    D3vector b( 0,18, 0);
    D3vector c( 0, 0,18);
    UnitCell uc(a,b,c);
    double ecut = 25.0;
    s.wf.resize(uc,uc,ecut);
    s.wf.set_nel(12*54);

    s.wf.randomize(1.e-4);
    s.wf.gram();
    cout << " ortho_error: " << s.wf.sd(0,0)->ortho_error() << endl;
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
