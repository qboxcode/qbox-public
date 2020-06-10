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
// testXMLGFPreprocessor.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
using namespace std;

#include "Context.h"
#include "Matrix.h"
#include "XMLGFPreprocessor.h"

int main(int argc, char** argv)
{
#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
  // extra scope to ensure that BlacsContext objects get destructed before
  // the MPI_Finalize call
  {
    if ( argc != 4 )
    {
      cout << "use: testXMLGFPreprocessor nprow npcol filename" << endl;
      return 1;
    }
    const int nr = atoi(argv[1]);
    const int nc = atoi(argv[2]);
    const char* const filename = argv[3];

    Context ctxt(MPI_COMM_WORLD,nr,nc); // context on which gfdata is defined
    DoubleMatrix gfdata(ctxt);
    string xmlcontent;

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int namelen;
    PMPI_Get_processor_name(processor_name,&namelen);
    cout << " Process " << gfdata.context().mype()
         << " on " << processor_name << endl;

    XMLGFPreprocessor xmlgfp;
    bool serial = true;
    xmlgfp.process(filename,gfdata,xmlcontent,serial);

#if 0
    // write all gfdata on file
    ofstream gf_file("gf.dat");
    gf_file << gfdata;
#endif

    if ( ctxt.onpe0() )
    {
      ofstream xmlfile("xml.dat");
      xmlfile << xmlcontent;
    }
  }
#if USE_MPI
  MPI_Finalize();
#endif
  return 0;
}
