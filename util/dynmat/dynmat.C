//
// dynmat.C: compute and diagonalize a dynamical matrix
// Forces are read from Qbox output
// The Qbox output should correspond to a sequence of calculations
// using symmetric finite displacements for all atoms in the x,y,z directions
//
// use: dynmat force.dat h Nat1 mass1 [Nat2 mass2] ...
// input_file: force.dat: forces from Qbox XML output file (collected with grep)
// h: displacement used in the force calculations (a.u.)
// Nat1: number of atoms of mass mass1
// Nat2: (optional) number of atoms of mass mass2
// (repeat the above for all atomic species)

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <valarray>
#include <vector>
using namespace std;

#include "Context.h"
#include "Matrix.h"

int main(int argc, char **argv)
{
  if ( argc < 5 )
  {
    cout << "use: dynmat force.dat h nat1 mass1 [nat2 mass2] ..." << endl;
    return 1;
  }
  char* infilename = argv[1];
  ifstream infile(infilename);

  const double h = atof(argv[2]);
  // cout << "h=" << h << endl;
  const double Ha2cm1 = 219474.65;

  vector<double> mass;
  for ( int iarg = 3; iarg < argc; iarg++, iarg++ )
  {
    // read number of atoms and mass for each species
    const int n = atoi(argv[iarg]);
    const double m = atof(argv[iarg+1]);
    // cout << n << " atoms of mass " << m << endl;
    for ( int i = 0; i < n; i++ )
      mass.push_back(1822.89*m);
  }
  const int nat = mass.size();

  // read forces
  double f[2][nat][3][nat][3];
  double fx,fy,fz;
  string dum;
  for ( int ia1 = 0; ia1 < nat; ia1++ )
  {
    for ( int idir1 = 0; idir1 < 3; idir1++ )
    {
      for ( int ishift = 0; ishift < 2; ishift++ )
      {
        for ( int ia2 = 0; ia2 < nat; ia2++ )
        {
          infile >> dum >> fx >> fy >> fz;
          while (infile.get() != '\n');
          f[ishift][ia1][idir1][ia2][0] = fx; 
          f[ishift][ia1][idir1][ia2][1] = fy; 
          f[ishift][ia1][idir1][ia2][2] = fz; 
        }
      }
    }
  }

  // compute matrix elements: centered finite difference
  // a[3*ia1+idir1][3*ia2+idir2] = 
  //  ( f[1][ia1][idir1][ia2][idir2] - f[0][ia1][idir1][ia2][idir2] ) / ( 2 h )

  Context ctxt;
  const int n = 3*nat;
  DoubleMatrix a(ctxt,n,n);

  for ( int ia1 = 0; ia1 < nat; ia1++ )
  {
    for ( int idir1 = 0; idir1 < 3; idir1++ )
    {
      for ( int ia2 = 0; ia2 < nat; ia2++ )
      {
        for ( int idir2 = 0; idir2 < 3; idir2++ )
        {
          int i = 3*ia1+idir1;
          int j = 3*ia2+idir2;
          double aij = ( f[1][ia1][idir1][ia2][idir2] - 
                         f[0][ia1][idir1][ia2][idir2] ) / ( 2.0 * h );
          aij = aij / sqrt(mass[ia1]*mass[ia2]);
          a[j*a.n()+i] = aij;
        }
      }
    }
  }

  // cout << a;
  valarray<double> w(a.m());
  DoubleMatrix asave(a);
  a.syevd('u',w);

  cout << " frequencies:" << endl;
  for ( int i = 0; i < a.n(); i++ )
  {
    if ( w[i] < 0.0 )
      cout << setw(8) << (int) (Ha2cm1 * sqrt(-w[i])) << " I" << endl;
    else
      cout << setw(8) << (int) (Ha2cm1 * sqrt(w[i])) << endl;
  }

  a = asave;
  a.syevd('l',w);
  cout << " frequencies:" << endl;
  for ( int i = 0; i < a.n(); i++ )
  {
    if ( w[i] < 0.0 )
      cout << setw(8) << (int) (Ha2cm1 * sqrt(-w[i])) << " I" << endl;
    else
      cout << setw(8) << (int) (Ha2cm1 * sqrt(w[i])) << endl;
  }
}
