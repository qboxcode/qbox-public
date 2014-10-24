//
// kpgen.C: generate a kpoint list for an arbitrary cell
// use: kpgen nx ny nz sx sy sz  a0x a0y a0z   a1x a1y a1z   a2x a2y a2z
//
// where a0x a0y a0z = first basis vector of the unit cell
//       a1x a1y a1z = second basis vector of the unit cell
//       a2x a2y a2z = third basis vector of the unit cell
//
//       nx,ny,nz: number of kpoints in each direction
//       sx,sy,sz: shift in each direction
//       shift: 0: symmetric set: boundary points not included
//              1: shifted set: boundary points included
//
#include<iostream>
#include<iomanip>
#include<vector>
#include<cassert>
#include<list>
#include "D3vector.h"
using namespace std;

int main(int argc, char** argv)
{
  if ( argc != 16 )
  {
    cerr << " use: " << argv[0] << " nx ny nz shiftx shifty shiftz {cell}"
         << endl;
    cerr << "      shift==0: symmetric set, zone boundary not included" << endl;
    cerr << "      shift==1: shifted set, zone boundary included" << endl;
    return 1;
  }
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);
  int sx = atoi(argv[4]);
  int sy = atoi(argv[5]);
  int sz = atoi(argv[6]);
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
    cerr << " shifts must be 0 or 1" << endl;
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
  for ( int i = nx-1; i >= 0; i-- )
  {
    for ( int j = ny-1; j >= 0; j-- )
    {
      for ( int k = nz-1; k >= 0; k-- )
      {
        kpint[0] = 2*i-nx+sx+1;
        kpint[1] = 2*j-ny+sy+1;
        kpint[2] = 2*k-nz+sz+1;
        kpint[3] = 1;
        kplist.push_back(kpint);
      }
    }
  }

  // fold kpoints back in the BZ
  for (list<vector<int> >::iterator i = kplist.begin();
       i != kplist.end(); i++)
  {
    kpint[0] = (*i)[0];
    kpint[1] = (*i)[1];
    kpint[2] = (*i)[2];
    kpint[3] = (*i)[3];
    bool done = false;
    while (!done)
    {
      done = true;
      D3vector kabs = kpint[0]/(2.0*nx) * b0 +
                      kpint[1]/(2.0*ny) * b1 +
                      kpint[2]/(2.0*nz) * b2;
      D3vector g;
      // check projection of kpoint along all 26 reciprocal lattice vectors
      // that are nearest g=0

      for ( int i0 = -1; i0 <= 1; i0++ )
      for ( int i1 = -1; i1 <= 1; i1++ )
      for ( int i2 = -1; i2 <= 1; i2++ )
      if ( !(i0 == 0 && i1 == 0 && i2 == 0) )
      {
        D3vector g = i0 * b0 + i1 * b1 + i2 * b2;
        if ( kabs*g >  (0.5 + 1.e-4) * g*g )
        {
          kpint[0]-= i0*2*nx;
          kpint[1]-= i1*2*ny;
          kpint[2]-= i2*2*nz;
          kabs = kpint[0]/(2.0*nx) * b0 +
                 kpint[1]/(2.0*ny) * b1 +
                 kpint[2]/(2.0*nz) * b2;
          done = false;
        }
      }
    } // while !done
    // kpint is now inside the BZ
    (*i)[0] = kpint[0];
    (*i)[1] = kpint[1];
    (*i)[2] = kpint[2];
  }

  // remove -k
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    // test if -k is in the set (and is not 0 0 0)
    kpint[0] = -(*i)[0];
    kpint[1] = -(*i)[1];
    kpint[2] = -(*i)[2];
    kpint[3] =  (*i)[3];
    if ( kpint[0]*kpint[0]+kpint[1]*kpint[1]+kpint[2]*kpint[2] != 0 )
    {
      // look for -k in the rest of the list
      list<vector<int> >::iterator j = find(i,kplist.end(),kpint);
      if ( j != kplist.end() )
      {
        kplist.erase(j);
        // double the weight of kpint
        (*i)[3] *= 2;
      }
    }
  }

  // check that sum of weights is one
  int sum = 0;
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    sum += (*i)[3];
  }
  assert(sum==nx*ny*nz);

  // output list
  // traverse list  backwards to have increasing indices
  // kpoints are output in reciprocal lattice coordinates
  cout.setf(ios::right,ios::adjustfield);
  cout << "# kpgen     " << nx << " " << ny << " " << nz << endl;
  cout << "# sx,sy,sz= " << sx << " " << sy << " " << sz << endl;
  cout << "# a0 = " << a0 << endl;
  cout << "# a1 = " << a1 << endl;
  cout << "# a2 = " << a2 << endl;
  cout << "# " << kplist.size() << " k-points" << endl;
  cout << " kpoint delete 0 0 0" << endl;
  for (list<vector<int> >::reverse_iterator i = kplist.rbegin();
       i != kplist.rend(); i++)
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout << " kpoint add "
         << setprecision(10)
         << setw(13) << (*i)[0]/(2.0*nx) << " "
         << setw(13) << (*i)[1]/(2.0*ny) << " "
         << setw(13) << (*i)[2]/(2.0*nz) << "   ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setprecision(14)
         << setw(16) << (*i)[3]/((double) nx*ny*nz)
         << endl;
  }

}
