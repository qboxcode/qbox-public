////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2010 The Regents of the University of California
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
// dynmat.cpp: compute and diagonalize a dynamical matrix
//
////////////////////////////////////////////////////////////////////////////////
// Forces are read from Qbox output
// The Qbox output should correspond to a sequence of calculations
// using symmetric finite displacements for all atoms in the x,y,z directions
// Displacements have an amplitude of +/- h
//
// use: dynmat force.dat h Nat1 mass1 [Nat2 mass2] ...
// input_file: force.dat: forces from Qbox XML output file (collected with grep)
// h: displacement used in the force calculations (a.u.)
// Nat1: number of atoms of mass mass1
// Nat2: (optional) number of atoms of mass mass2
// (repeat the above for all atomic species)
// Species masses must be given in the order in which they appear in the
// <atomset> element in Qbox output

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

extern "C"
void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
  const int *lda, double *w, double *wrk, int *lwrk, int *info);

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

  const int n = 3*nat;
  valarray<double> a(n*n);

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
          a[j*n+i] = aij;
        }
      }
    }
  }

  // cout << a;
  valarray<double> w_upper(n),w_lower(n),w_sym(n);
  valarray<double> asave(a);
  char jobz = 'n';
  char uplo = 'u';
  int lwrk = 3*n;
  valarray<double> wrk(lwrk);

  // use the upper part of the dynamical matrix
  int info;
  dsyev_(&jobz,&uplo,&n,&a[0],&n,&w_upper[0],&wrk[0],&lwrk,&info);
  assert(info==0);

  a = asave;
  // use the lower part of the dynamical matrix
  uplo = 'l';
  dsyev_(&jobz,&uplo,&n,&a[0],&n,&w_lower[0],&wrk[0],&lwrk,&info);
  assert(info==0);

  a = asave;
  // symmetrize the matrix a
  for ( int i = 0; i < n; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      double aij = 0.5 * ( a[j*n+i] + a[i*n+j] );
      a[j*n+i] = aij;
      a[i*n+j] = aij;
    }
  }

  // compute eigenvectors
  jobz = 'v';
  dsyev_(&jobz,&uplo,&n,&a[0],&n,&w_sym[0],&wrk[0],&lwrk,&info);
  assert(info==0);

  cout << " frequencies (upper/lower/sym):" << endl;
  for ( int i = 0; i < n; i++ )
  {
    if ( w_upper[i] < 0.0 )
      cout << setw(8) << (int) (Ha2cm1 * sqrt(-w_upper[i])) << " I";
    else
      cout << setw(8) << (int) (Ha2cm1 * sqrt(w_upper[i]));

    if ( w_lower[i] < 0.0 )
      cout << setw(8) << (int) (Ha2cm1 * sqrt(-w_lower[i])) << " I";
    else
      cout << setw(8) << (int) (Ha2cm1 * sqrt(w_lower[i]));

    if ( w_sym[i] < 0.0 )
      cout << setw(8) << (int) (Ha2cm1 * sqrt(-w_sym[i])) << " I";
    else
      cout << setw(8) << (int) (Ha2cm1 * sqrt(w_sym[i]));
    cout << endl;
  }

  ofstream vecfile("dynmat_eigvec.dat");
  for ( int j = 0; j < n; j++ )
  {
    vecfile << "# mode " << j+1 << "  frequency = ";
    if ( w_sym[j] < 0.0 )
      vecfile << setw(8) << (int) (Ha2cm1 * sqrt(-w_sym[j])) << " I";
    else
      vecfile << setw(8) << (int) (Ha2cm1 * sqrt(w_sym[j]));
    vecfile << " cm-1" << endl;
    vecfile << setprecision(8);
    vecfile.setf(ios::fixed, ios::floatfield);
    vecfile.setf(ios::right, ios::adjustfield);
    for ( int i = 0; i < n; i+=3 )

      vecfile << i/3+1 << " "
              << setw(12) << a[j*n+i] << " "
              << setw(12) << a[j*n+i+1] << " "
              << setw(12) << a[j*n+i+2] << endl;
  }
}
