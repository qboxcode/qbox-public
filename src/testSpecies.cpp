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
// testSpecies.cpp
//
// use: testSpecies uri
//

#include "Species.h"
#include "SpeciesReader.h"
#include <iostream>
#include <cassert>
#include <string>
using namespace std;

#include <mpi.h>
int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  {

  if ( argc != 2 )
  {
    cerr << "use: testSpecies uri" << endl;
    return 1;
  }

  Species s("unknown_name");

  SpeciesReader rdr;

  string uri(argv[1]);

  cout << " s.uri() = " << s.uri() << endl;
  cout << " testSpecies: invoking SpeciesReader::uri_to_species:" << endl;
  rdr.uri_to_species(uri,s);
  cout << " SpeciesReader::uri_to_species done" << endl;

  const double rcps = 1.0;

  try
  {
    s.initialize(rcps);
  }
  catch ( SpeciesInitException& e )
  {
    cerr << " Exception in Species initialization: " << e.msg << endl;
    throw;
  }
  cout << s;

  cout << " testSpecies: output of species done" << endl;

  cout << " testSpecies: invoking s.info(cout):" << endl;
  s.info(cout);

  cout << " testSpecies: testing SpeciesReader::uri_to_string: " << endl;
  string xmlstr;
  rdr.uri_to_string(uri,"unknown",xmlstr);
  cout << xmlstr;

#if 0

  if ( ctxt.onpe0() )
  {

  double dr = 0.01;
  double dg = 0.02;
  double v,dv;

  int n = 1500;

  for ( int l = 0; l <= s.lmax(); l++ )
  {
    cout << "# " << n << " Vps(l=" << l << ",r) " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.vpsr(l,r,v);
      cout << r << " " << v << endl;
    }
    cout << endl << endl;

    cout << "# " << n << " phi(l=" << l << ",r) " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      double val;
      s.phi(l,r,val);
      cout << r << " " << val << endl;
    }
    cout << endl << endl;

    cout << "# " << n << " dVps(l=" << l << ",r)/dr " << endl;
    for ( int i = 0; i < n; i++ )
    {
      double r = i * dr;
      s.dvpsr(l,r,v,dv);
      cout << r << " " << dv << endl;
    }
    cout << endl << endl;
  }

  cout << "# " << n << " Vloc(g) " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.vlocg(g,v);
    cout << g << " " << v << endl;
  }
  cout << endl << endl;

  cout << "# " << n << " dVloc(g)/dg " << endl;
  for ( int i = 0; i < n; i++ )
  {
    double g = i * dg;
    s.dvlocg(g,v,dv);
    cout << g << " " << dv << endl;
  }
  cout << endl << endl;

  for ( int l = 0; l <= s.lmax(); l++ )
  {
    if ( l != s.llocal() )
    {
      cout << "# " << n << " Vnl(l=" << l << ",g) " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.vnlg(l,g,v);
        cout << g << " " << v << endl;
      }
      cout << endl << endl;

      cout << "# " << n << " dVnl(l=" << l << ",g)/dg " << endl;
      for ( int i = 0; i < n; i++ )
      {
        double g = i * dg;
        s.dvnlg(l,g,v,dv);
        cout << g << " " << dv << endl;
      }
      cout << endl << endl;
    }
  }

  } // onpe0
#endif

  }
  MPI_Finalize();
  return 0;
}
