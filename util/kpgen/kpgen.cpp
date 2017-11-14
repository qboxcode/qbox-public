////////////////////////////////////////////////////////////////////////////////
//
// kpgen.cpp: generate a kpoint list for an arbitrary cell
// use: kpgen nx ny nz sx sy sz  a0x a0y a0z   a1x a1y a1z   a2x a2y a2z
//
// where a0x a0y a0z = first basis vector of the unit cell
//       a1x a1y a1z = second basis vector of the unit cell
//       a2x a2y a2z = third basis vector of the unit cell
//
//       nx,ny,nz: number of kpoints in each direction
//       sx,sy,sz: shift in each direction (floating point in [0,1) )
//
// note: from version 2.1 on, sx=sy=sz=0 generates a symmetric set
// A symmetric set includes the Gamma (k=0) point
//
// note: using non-zero shifts with non-rectangular cells breaks symmetry
// non-zero shifts should only be used with rectangular cells.
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cstdlib>
#include "D3vector.h"
using namespace std;

const char *const version = "2.1";

// largest shift when scanning reciprocal space
const int shmax=2;

////////////////////////////////////////////////////////////////////////////////
// comparison function for k-points
// two k-points are equal if k1 = k2+G for any G in the first shells
bool equals(D3vector k1, D3vector k2, D3vector b0, D3vector b1, D3vector b2)
{
  for ( int i0 = -shmax; i0 <= shmax; i0++ )
    for ( int i1 = -shmax; i1 <= shmax; i1++ )
      for ( int i2 = -shmax; i2 <= shmax; i2++ )
      {
        if ( length(k1-k2-i0*b0-i1*b1-i2*b2) < 1.e-5 )
          return true;
      }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// check if a k-point is in the BZ
bool in_BZ(D3vector k, D3vector b0, D3vector b1, D3vector b2)
{
  // check projection of kpoint k along all reciprocal lattice vectors
  // in the first shmax shells
  // use a shift in an irrational direction epsilon*(1,M_PI,M_LN2)
  // to avoid including zone boundary equivalent vectors
  const double epsilon = 1.e-6;
  D3vector kshifted = k + epsilon * D3vector(1.0,M_PI,M_LN2);
  for ( int i0 = -shmax; i0 <= shmax; i0++ )
    for ( int i1 = -shmax; i1 <= shmax; i1++ )
      for ( int i2 = -shmax; i2 <= shmax; i2++ )
        if ( !((i0 == 0) && (i1 == 0) && (i2 == 0)) )
        {
          D3vector g = i0 * b0 + i1 * b1 + i2 * b2;
          if ( kshifted*g >  0.5 * g*g )
            return false;
        }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  cout << "# kpgen " << version << endl;
  if ( argc != 16 )
  {
    cerr << " use: " << argv[0] << " nx ny nz shiftx shifty shiftz {cell}"
         << endl;
    return 1;
  }
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);
  double sx = atof(argv[4]);
  double sy = atof(argv[5]);
  double sz = atof(argv[6]);
  if ( nx <= 0 || ny <= 0 || nz <=0 )
  {
    cerr << " use: " << argv[0] << " nx ny nz shiftx shifty shiftz {cell}"
         << endl;
    cerr << " nx, ny, nz must be positive" << endl;
    return 1;
  }
  if ( sx < 0 || sx > 1 ||
       sy < 0 || sy > 1 ||
       sz < 0 || sz > 1 )
  {
    cerr << " use: " << argv[0] << " nx ny nz shiftx shifty shiftz {cell}"
         << endl;
    cerr << " shifts must be in [0,1]" << endl;
    return 1;
  }
  D3vector a0(atof(argv[7]),atof(argv[8]),atof(argv[9]));
  D3vector a1(atof(argv[10]),atof(argv[11]),atof(argv[12]));
  D3vector a2(atof(argv[13]),atof(argv[14]),atof(argv[15]));
  const double volume = a0 * ( a1 ^ a2 );
  if ( volume == 0.0 )
  {
    cout << " cell volume is zero" << endl;
    return 1;
  }
  double fac = 2.0 * M_PI / volume;
  D3vector b0 = fac * a1 ^ a2;
  D3vector b1 = fac * a2 ^ a0;
  D3vector b2 = fac * a0 ^ a1;

  // check if shifts are used with a non-rectangular cell
  if ( sx > 0 || sy > 0 || sz > 0 )
  {
    if ( (fabs(a0*a1) > 1.e-6) ||
         (fabs(a1*a2) > 1.e-6) ||
         (fabs(a0*a2) > 1.e-6) )
    {
      cout << " warning: non-zero shifts with non-rectangular cell"
           << " may break symmetry" << endl;
    }
  }

  vector<D3vector> kp;
  vector<D3vector> kpfrac;
  vector<double> weight;

  // scan volume enclosing the BZ
  for ( int ii = -shmax; ii <= shmax; ii++ )
  for ( int jj = -shmax; jj <= shmax; jj++ )
  for ( int kk = -shmax; kk <= shmax; kk++ )
  for ( int i = 0; i < nx; i++ )
  {
    for ( int j = 0; j < ny; j++ )
    {
      for ( int k = 0; k < nz; k++ )
      {
        int kpint0 = ii*2*nx + 2*i-nx+1;
        int kpint1 = jj*2*ny + 2*j-ny+1;
        int kpint2 = kk*2*nz + 2*k-nz+1;

        double kv0 = (kpint0 + sx + (nx%2+1))/(2.0*nx);
        double kv1 = (kpint1 + sy + (ny%2+1))/(2.0*ny);
        double kv2 = (kpint2 + sz + (nz%2+1))/(2.0*nz);

        D3vector kv = kv0*b0 + kv1*b1 + kv2*b2;

        if ( in_BZ(kv,b0,b1,b2) )
        {
          kp.push_back(kv);
          kpfrac.push_back(D3vector(kv0,kv1,kv2));
          weight.push_back(1.0);
        }
      }
    }
  }
  int total_weight = kp.size();

  // check for equivalent vectors
  // count vectors that are equivalent to k+G
  cout.setf(ios::fixed,ios::floatfield);
#if 0
  int nequiv = 0;
  for ( int i = 0; i < kp.size(); i++ )
    for ( int j = i+1; j < kp.size(); j++ )
      if ( equals(kp[i],kp[j],b0,b1,b2) )
      {
        nequiv++;
        //cout << setprecision(3)
        //     << kpfrac[i] << " equiv " << kpfrac[j] << endl;

      }
  if ( nequiv != 0 )
  {
    // there should not be any equivalent points as k=k+G
    cerr << nequiv << " error: equivalent points (k=k+G)" << endl;
    //return 1;
  }
#endif

  // count vectors that are equivalent to -k+G
#if 0
  int nequivm = 0;
  for ( int i = 0; i < kp.size(); i++ )
    for ( int j = i+1; j < kp.size(); j++ )
      if ( equals(kp[i],-kp[j],b0,b1,b2) )
      {
        nequivm++;
        //cout << setprecision(3)
        //     << kpfrac[i] << " equiv " << kpfrac[j] << endl;
      }
  cout << nequivm << " equivalent points (k=-k+G)" << endl;
#endif

  // check duplicates
#if 0
  int ndup = 0;
  for ( int i = 0; i < kp.size(); i++ )
    for ( int j = i+1; j < kp.size(); j++ )
      if ( length(kp[i]-kp[j]) < 1.e-5 )
      {
        ndup++;
        cout << setprecision(3)
             << kpfrac[i] << " duplicate of " << kpfrac[j] << endl;
      }
  cout << "# " << ndup << " duplicates" << endl;
#endif

#if 1
  // reassign weight from (k,-k+G) equivalent points
  for ( int i = 0; i < kp.size(); i++ )
  {
    if ( weight[i] != 0.0 )
    {
      for ( int j = i+1; j < kp.size(); j++ )
      {
        if ( equals(kp[i],-kp[j],b0,b1,b2) )
        {
          // reassign the weight of kp[j] to kp[i]
          weight[i] += weight[j];
          weight[j] = 0.0;
        }
      }
    }
  }
#endif

  // count k points with non-zero weight
  int nkp = 0;
  for ( int i = 0; i < weight.size(); i++ )
    if ( weight[i] > 0.0 ) nkp++;

  // output list
  // kpoints are output in reciprocal lattice coordinates
  cout.setf(ios::right,ios::adjustfield);
  cout << "# nx,ny,nz: " << nx << " " << ny << " " << nz << endl;
  cout << "# sx,sy,sz: " << sx << " " << sy << " " << sz << endl;
  cout << "# a0: " << a0 << endl;
  cout << "# a1: " << a1 << endl;
  cout << "# a2: " << a2 << endl;
  cout << "# b0/(2pi): " << b0/(2*M_PI) << endl;
  cout << "# b1/(2pi): " << b1/(2*M_PI) << endl;
  cout << "# b2/(2pi): " << b2/(2*M_PI) << endl;

  cout << "# " << nkp << " k-points" << endl;
  cout << " kpoint delete 0 0 0" << endl;

  // print list backward to have increasing x components
  for ( int i = kpfrac.size()-1; i >= 0; i-- )
  {
    if ( weight[i] != 0.0 )
    {
      cout.setf(ios::fixed,ios::floatfield);
      // compute projections along b0, b1, b2
      double kx = kpfrac[i].x;
      double ky = kpfrac[i].y;
      double kz = kpfrac[i].z;
      // print -k to have positive coefficients in output
      // change sign only if component is non-zero to avoid -0.00 in output
#if 1
      if ( kx == 0.0 )
      {
        if ( ky == 0.0 )
        {
          if ( kz < 0.0 ) kz = -kz;
        }
        else if ( ky < 0.0 )
        {
          ky = -ky;
          if ( kz != 0 ) kz = -kz;
        }
      }
      else if ( kx < 0.0 )
      {
        kx = -kx;
        if ( ky != 0 ) ky = -ky;
        if ( kz != 0 ) kz = -kz;
      }
#endif
      double w = weight[i]/((double) total_weight);

      cout << " kpoint add "
           << setprecision(10)
           << setw(13) << kx << " "
           << setw(13) << ky << " "
           << setw(13) << kz << "   ";
      cout.setf(ios::scientific,ios::floatfield);
      cout << setprecision(14)
           << setw(16) << w << endl;
    }
  }
}
