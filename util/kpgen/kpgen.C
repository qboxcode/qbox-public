//
// kpgen.C: generate a kpoint list for an arbitrary cell
// use: kpgen nx ny nz sx sy sz  a0x a0y a0z   a1x a1y a1z   a2x a2y a2z
//
// where a0x a0y a0z = first basis vector of the unit cell
//       a1x a1y a1z = second basis vector of the unit cell
//       a2x a2y a2z = third basis vector of the unit cell
//
//       nx,ny,nz: number of kpoints in each direction
//       sx,sy,sz: shift in each direction (floating point)
//       shift: 0: symmetric set: boundary points not included
// even-numbered sets do not include the gamma point
// odd-numbered sets include the gamma point
//
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cassert>
#include<cstdlib>
#include<list>
#include "D3vector.h"
using namespace std;

// ib_BZ: test if the vector k is in the BZ defined by b0,b1,b2
bool in_BZ(D3vector k, D3vector b0, D3vector b1, D3vector b2)
{
  const double epsilon = 1.e-8;
  D3vector g;
  // check projection of kpoint k along all 26 reciprocal lattice vectors
  // that are nearest g=0
  // use a shift by epsilon*(1,1,1) to avoid including zone boundary
  // equivalent vectors

  bool in_bz = true;
  D3vector kshifted = k + epsilon * D3vector(1.0,1.0,1.0);
  for ( int i0 = -1; i0 <= 1; i0++ )
    for ( int i1 = -1; i1 <= 1; i1++ )
      for ( int i2 = -1; i2 <= 1; i2++ )
      if ( !(i0 == 0 && i1 == 0 && i2 == 0) )
      {
        D3vector g = i0 * b0 + i1 * b1 + i2 * b2;
        if ( kshifted*g >  0.5 * g*g )
          in_bz = false;
      }
  return in_bz;
}

int main(int argc, char** argv)
{
  cout << "# kpgen " << endl;
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

  list<vector<int> > kplist;
  vector<int> kpint(4);

  // scan volume enclosing the BZ
  for ( int ii = -2; ii <= 2; ii++ )
  for ( int jj = -2; jj <= 2; jj++ )
  for ( int kk = -2; kk <= 2; kk++ )
  for ( int i = 0; i < nx; i++ )
  {
    for ( int j = 0; j < ny; j++ )
    {
      for ( int k = 0; k < nz; k++ )
      {
        kpint[0] = ii*2*nx + 2*i-nx+1;
        kpint[1] = jj*2*ny + 2*j-ny+1;
        kpint[2] = kk*2*nz + 2*k-nz+1;
        kpint[3] = 1;

        D3vector kv = ( (kpint[0] + sx)/(2.0*nx) ) * b0 +
                      ( (kpint[1] + sy)/(2.0*ny) ) * b1 +
                      ( (kpint[2] + sz)/(2.0*nz) ) * b2;

        if ( in_BZ(kv,b0,b1,b2) )
          kplist.push_back(kpint);
      }
    }
  }

  int total_weight = kplist.size();

#if 1
  // remove -k
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    // test if -k is in the set (and is not 0 0 0)
    kpint[0] = (*i)[0];
    kpint[1] = (*i)[1];
    kpint[2] = (*i)[2];
    kpint[3] = (*i)[3];
    D3vector ki = (((*i)[0]+sx)/(2.0*nx)) * b0 +
                  (((*i)[1]+sy)/(2.0*ny)) * b1 +
                  (((*i)[2]+sz)/(2.0*nz)) * b2;
    if ( length(ki) != 0.0 )
    {
      // look for -k in the rest of the list
      for ( list<vector<int> >::iterator j = i; j != kplist.end(); j++ )
      {
        D3vector kj = (((*j)[0]+sx)/(2.0*nx)) * b0 +
                      (((*j)[1]+sy)/(2.0*ny)) * b1 +
                      (((*j)[2]+sz)/(2.0*nz)) * b2;
        if ( length(ki+kj) < 1.e-5 )
        {
          // transfer weight to (*i)
          (*i)[3] += (*j)[3];
          (*j)[3] = 0;
        }
      }
    }
  }
#endif

#if 1
  // remove duplicate points
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    D3vector ki = (((*i)[0]+sx)/(2.0*nx)) * b0 +
                  (((*i)[1]+sy)/(2.0*ny)) * b1 +
                  (((*i)[2]+sz)/(2.0*nz)) * b2;
    // look for duplicate points in the rest of the list
    for ( list<vector<int> >::iterator j = i; j != kplist.end(); j++ )
    {
      if ( j != i )
      {
        D3vector kj = (((*j)[0]+sx)/(2.0*nx)) * b0 +
                      (((*j)[1]+sy)/(2.0*ny)) * b1 +
                      (((*j)[2]+sz)/(2.0*nz)) * b2;
        if ( length(ki-kj) < 1.e-5 )
        {
          // transfer the weight of kj to ki
          (*i)[3] += (*j)[3];
          (*j)[3] = 0;
        }
      }
    }
  }
#endif


#if 1
  // check that the sum of weights is one
  int sum = 0;
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    sum += (*i)[3];
  }
  assert(sum==total_weight);
#endif

#if 1
  // remove elements with zero weight
  // For a description of safe erase, see S. Meyers, Effective STL, item 9
  for (list<vector<int> >::iterator i = kplist.begin();
       i != kplist.end(); /* nothing */)
  {
    int w = (*i)[3];
    if ( w == 0 )
      kplist.erase(i++);
    else
      ++i;
  }
#endif

#if 1
  // output list
  // change the sign of the k vector if the first element is negative
  // traverse list  backwards to have increasing indices
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

  cout << "# " << kplist.size() << " k-points" << endl;
  cout << " kpoint delete 0 0 0" << endl;
  for (list<vector<int> >::reverse_iterator i = kplist.rbegin();
       i != kplist.rend(); i++)
  {
    D3vector ki = ( ((*i)[0]+sx)/(2.0*nx) ) * b0 +
                  ( ((*i)[1]+sy)/(2.0*ny) ) * b1 +
                  ( ((*i)[2]+sz)/(2.0*nz) ) * b2;
    cout.setf(ios::fixed,ios::floatfield);
    double kx = ((*i)[0]+sx)/(2.0*nx);
    double ky = ((*i)[1]+sy)/(2.0*ny);
    double kz = ((*i)[2]+sz)/(2.0*nz);
    double w = (*i)[3]/((double) total_weight);
    if ( ki.x < 0.0 )
    {
      kx = -kx;
      // next lines: test before changing sign to avoid -0.0
      if ( ky != 0.0 ) ky = -ky;
      if ( kz != 0.0 ) kz = -kz;
    }

    cout << " kpoint add "
         << setprecision(10)
         << setw(13) << kx << " "
         << setw(13) << ky << " "
         << setw(13) << kz << "   ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setprecision(14)
         << setw(16) << w << endl;
  }
#endif

#if 1
  // test the k-point set
  // compute the numerical integral of the function exp(ikR) for
  // a set of vectors R. The integral should be zero, except for R=0 and
  // vectors R of large norm (larger than the nx,ny,nz parameters used).
  double minlen = 1.e10;
  for ( int ii = -nx; ii <= nx; ii++ )
  {
    for ( int jj = -ny; jj <= ny; jj++ )
    {
      for ( int kk = -nz; kk <= nz; kk++ )
      {
        D3vector R = ii * a0 + jj * a1 + kk * a2;
        double len = length(R);
        double sum = 0.0;
        for (list<vector<int> >::iterator ikp = kplist.begin();
             ikp != kplist.end(); ikp++)
        {
          D3vector k = ( ((*ikp)[0]+sx)/(2.0*nx) ) * b0 +
                       ( ((*ikp)[1]+sy)/(2.0*ny) ) * b1 +
                       ( ((*ikp)[2]+sz)/(2.0*nz) ) * b2;
          double w = (*ikp)[3]/((double) total_weight);
          sum += w * cos(k*R);
        }

        if ( len != 0 && fabs(sum) > 1.e-6 )
          minlen = min(minlen,len);
      }
    }
  }
  cerr << " smallest R with non-zero sum has length " << minlen << endl;
#endif
}


