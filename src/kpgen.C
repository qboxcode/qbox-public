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

  bool in_bz = true;
  for ( int i0 = -1; i0 <= 1; i0++ )
    for ( int i1 = -1; i1 <= 1; i1++ )
      for ( int i2 = -1; i2 <= 1; i2++ )
      if ( !(i0 == 0 && i1 == 0 && i2 == 0) )
      {
        D3vector g = i0 * b0 + i1 * b1 + i2 * b2;
        if ( k*g >  0.5 * g*g + epsilon )
          in_bz = false;
      }
  return in_bz;
}

int main(int argc, char** argv)
{
  cout << "# kpgen-1.0" << endl;
  if ( argc != 16 )
  {
    cerr << " use: " << argv[0] << " nx ny nz shiftx shifty shiftz {cell}"
         << endl;
    cerr << "      shift==0: symmetric set, zone boundary not included" << endl;
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
  for ( int kk = -1; kk <= 2; kk++ )
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

        D3vector k = kpint[0]/(2.0*nx) * b0 +
                     kpint[1]/(2.0*ny) * b1 +
                     kpint[2]/(2.0*nz) * b2;

        if ( in_BZ(k,b0,b1,b2) )
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
    if ( kpint[0]*kpint[0]+kpint[1]*kpint[1]+kpint[2]*kpint[2] != 0 )
    {
      // look for -k in the rest of the list
      for ( list<vector<int> >::iterator j = i; j != kplist.end(); j++ )
      {
        if ( (*j)[0]==-kpint[0] && (*j)[1]==-kpint[1] && (*j)[2]==-kpint[2] )
        {
          // transfer weight to (*i)
          (*i)[3] += (*j)[3];
          (*j)[3] = 0;
          //cout << " erasing  " << "(" << kpint[0] << ","
          //     << kpint[1] << "," << kpint[2] << ") == -("
          //     << (*j)[0] << "," << (*j)[1] << "," << (*j)[2] << ")" << endl;
        }
      }
    }
  }
#endif

#if 1
  // remove equivalent points
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    D3vector ki = (*i)[0]/(2.0*nx) * b0 +
                  (*i)[1]/(2.0*ny) * b1 +
                  (*i)[2]/(2.0*nz) * b2;
    // look for equivalent points in the rest of the list
    for ( list<vector<int> >::iterator j = i; j != kplist.end(); j++ )
    {
      if ( j != i )
      {
        D3vector kj = (*j)[0]/(2.0*nx) * b0 +
                      (*j)[1]/(2.0*ny) * b1 +
                      (*j)[2]/(2.0*nz) * b2;
        if ( length(ki-kj) < 1.e-5 )
        {
          // transfer the weight of kj to ki
          (*i)[3] += (*j)[3];
          (*j)[3] = 0;
          //cout << " erasing equivalent point " << "(" << (*j)[0] << ","
          //     << (*j)[1] << "," << (*j)[2] << ") == ("
          //     << (*i)[0] << "," << (*i)[1] << "," << (*i)[2] << ")" << endl;
        }
      }
    }
  }
#endif

  // remove elements with zero weight
  for (list<vector<int> >::iterator i = kplist.begin();
       i != kplist.end(); /* nothing */ )
  {
    if ( (*i)[3] == 0 )
    {
      kplist.erase(i++);
    }
    else
    {
      i++;
    }
  }

#if 1
  // make first index positive
  for (list<vector<int> >::iterator i = kplist.begin(); i != kplist.end(); i++)
  {
    D3vector ki = (*i)[0]/(2.0*nx) * b0 +
                  (*i)[1]/(2.0*ny) * b1 +
                  (*i)[2]/(2.0*nz) * b2;
    if ( ki.x < 0 )
    {
      (*i)[0] *= -1;
      (*i)[1] *= -1;
      (*i)[2] *= -1;
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
  // output list
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
    cout.setf(ios::fixed,ios::floatfield);
    cout << " kpoint add "
         << setprecision(10)
         << setw(13) << ((*i)[0]+sx)/(2.0*nx) << " "
         << setw(13) << ((*i)[1]+sy)/(2.0*ny) << " "
         << setw(13) << ((*i)[2]+sz)/(2.0*nz) << "   ";
    cout.setf(ios::scientific,ios::floatfield);
    cout << setprecision(14)
         << setw(16) << (*i)[3]/((double) total_weight)
         << endl;
  }

  // output list in absolute coordinates for plot
  ofstream plotfile("kpoint.plt");
  for (list<vector<int> >::reverse_iterator i = kplist.rbegin();
       i != kplist.rend(); i++)
  {
    D3vector k = (((*i)[0]+sx)/(2.0*nx)) * b0 / (2*M_PI) +
                 (((*i)[1]+sy)/(2.0*ny)) * b1 / (2*M_PI) +
                 (((*i)[2]+sz)/(2.0*nz)) * b2 / (2*M_PI);
    plotfile.setf(ios::fixed,ios::floatfield);
    plotfile << setprecision(8)
         << setw(13) << k.x << " "
         << setw(13) << k.y << " "
         << setw(13) << k.z << endl;

  }
#endif
}
