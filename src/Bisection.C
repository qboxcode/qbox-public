////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// Bisection.C
//
////////////////////////////////////////////////////////////////////////////////
#include "Bisection.h"
#include "VectorLess.h"
#include <bitset>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
int walsh(int l, int n, int i)
{
  // walsh function at level l at position i in a grid of n points
  assert(i>=0);
  assert(i<n);
  if ( l == 0 )
  {
    assert(n%2==0);
    if ( i >= n/2 ) return 0;
    else return 1;
  }
  else if ( l == 1 )
  {
    assert(n%4==0);
    if ( i >= n/4 && i < (3*n)/4 ) return 0;
    else return 1;
  }
  else if ( l == 2 )
  {
    assert(n%8==0);
    if ( (i >= n/8 && i < 3*n/8) || (i >= 5*n/8 && i < 7*n/8) ) return 0;
    else return 1;
  }
  else if ( l == 3 )
  {
    assert(n%16==0);
    if ( (i >= n/16 && i < 3*n/16) ||
         (i >= 5*n/16 && i < 7*n/16) ||
         (i >= 9*n/16 && i < 11*n/16) ||
         (i >= 13*n/16 && i < 15*n/16) ) return 0;
    else return 1;
  }
  else if ( l == 4 )
  {
    assert(n%32==0);
    if ( (i >= n/32 && i < 3*n/32) ||
         (i >= 5*n/32 && i < 7*n/32) ||
         (i >= 9*n/32 && i < 11*n/32) ||
         (i >= 13*n/32 && i < 15*n/32) ||
         (i >= 17*n/32 && i < 19*n/32) ||
         (i >= 21*n/32 && i < 23*n/32) ||
         (i >= 25*n/32 && i < 27*n/32) ||
         (i >= 29*n/32 && i < 31*n/32) ) return 0;
    else return 1;
  }
  else
    assert(false);
}

////////////////////////////////////////////////////////////////////////////////
Bisection::Bisection(Sample& s, int nlevels[3]) : s_(s)
{
  // localization vectors are long int
  assert(sizeof(long int) >= 4);

  Wavefunction& wf = s_.wf;
  gcontext_ = wf.sd(0,0)->context();
  // assume that nst is the same for both spins
  // (some arrays below are allocated with same size for both spins)
  if ( wf.nspin() > 1 )
    assert(wf.sd(0,0)->nst() == wf.sd(1,0)->nst());

  // number of bisection levels in each direction
  nlevels_[0]=nlevels[0];
  nlevels_[1]=nlevels[1];
  nlevels_[2]=nlevels[2];

  // number of subdivisions required in each direction
  // ndiv = 2^nlevel
  ndiv_[0] = 1 << nlevels[0];
  ndiv_[1] = 1 << nlevels[1];
  ndiv_[2] = 1 << nlevels[2];

  // largest number of levels
  nlevelsmax_=max(nlevels[0],max(nlevels[1],nlevels[2]));

  // real-space grid size for wave functions
  np_[0] = wf.sd(0,0)->basis().np(0);
  np_[1] = wf.sd(0,0)->basis().np(1);
  np_[2] = wf.sd(0,0)->basis().np(2);
  // adapt the grid dimensions to the levels of bisection
  for ( int idim = 0; idim < 3; idim++ )
  {
    for ( int j=1; j<=nlevels[idim]; j++ )
    {
      int base = 1 << j;
      if ( np_[idim] % base != 0 ) np_[idim] += base/2;
    }
  }

  // number of grid points of augmented grid for normalization
  ft_ = new FourierTransform(wf.sd(0,0)->basis(),np_[0],np_[1],np_[2]);
  np01_ = np_[0]*np_[1];
  np2loc_ = ft_->np2_loc();
  np012loc_ = ft_->np012loc();

  // xy projector index: xy_proj[i+j*np0] is the xy projector
  // associated with a point located at i + j*np0 + k*np0*np1
  xy_proj_.resize(np01_);

  for ( int i = 0; i < np_[0]; i++ )
    for ( int j = 0; j < np_[1]; j++ )
    {
      int i_slice_x= ( i * ndiv_[0] ) / np_[0];
      int i_slice_y= ( j * ndiv_[1] ) / np_[1];
      int xy_proj_index = i_slice_x + ndiv_[0] * i_slice_y;
      xy_proj_[i+np_[0]*j] = xy_proj_index;
    }

  // nmat_: number of A matrices
  nmat_ = nlevels[0] + nlevels[1] + nlevels[2];

  // allocate A matrices
  amat_.resize(nmat_);
  const ComplexMatrix &c = wf.sd(0,0)->c();
  for ( int i = 0; i < nmat_; i++ )
    amat_[i] = new DoubleMatrix(c.context(),c.n(),c.n(),c.nb(),c.nb());
  // allocate rotation matrix
  u_ = new DoubleMatrix(c.context(),c.n(),c.n(),c.nb(),c.nb());

  // matrices of real space wave functions in subdomains
  rmat_.resize( ndiv_[0]*ndiv_[1] );
  {
    // xyproj_rsize: real-space size of xy projectors
    vector<int> xyproj_rsize( ndiv_[0]*ndiv_[1] , 0 );
    for ( int ixy = 0; ixy < np01_; ixy++ )
      xyproj_rsize[xy_proj_[ixy]] += np2loc_;

    // max_xyproj_rsize: maximum real-space size of xy projectors
    vector<int> max_xyproj_rsize( ndiv_[0]*ndiv_[1] , 0 );
    MPI_Allreduce((void*)&xyproj_rsize[0] , (void*)&max_xyproj_rsize[0],
      (int)xyproj_rsize.size(), MPI_INT ,MPI_MAX , gcontext_.comm());

    // allocate matrices rmat_[i]
    for ( int i = 0; i < rmat_.size(); i++ )
    {
      int n = wf.sd(0,0)->c().n();
      int nb = wf.sd(0,0)->c().nb();
      int m = max_xyproj_rsize[i] * wf.sd(0,0)->c().context().nprow();
      int mb = max_xyproj_rsize[i];
      rmat_[i] = new DoubleMatrix(wf.sd(0,0)->c().context(),m,n,mb,nb);
      rmat_[i]->clear();
    }
  }

  localization_.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    localization_[ispin].resize(wf.sd(ispin,0)->nst());
  }
}

////////////////////////////////////////////////////////////////////////////////
Bisection::~Bisection(void)
{
  delete ft_;
  for ( int i = 0; i < nmat_; i++ ) delete amat_[i];
  delete u_;
  for ( int i = 0; i < rmat_.size(); i++ ) delete rmat_[i];
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::localize(Wavefunction& wf, double epsilon)
{
  map<string,Timer> tmap;

  // assert that we have only gamma point wave functions as jade is
  // not implemented for complex matrices
  assert( wf.nkp()==1 && wf.sd(0,0)->basis().real() );

  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    ComplexMatrix& c  = wf.sd(ispin,0)->c();
    const int nst = wf.sd(ispin,0)->nst();
    const int nstloc = wf.sd(ispin,0)->nstloc();
    const vector<double>& occ = wf.sd(ispin,0)->occ();
    assert(occ.size() == nst);

#ifdef TIMING
    tmap["wf FFT"].start();
#endif
    // Fourier transform states and save real-space values in
    // rmat_ matrices
    vector<complex<double> > wftmp(np012loc_,0.0);
    for ( int n = 0; n < nstloc; n++ )
    {
      ft_->backward(c.cvalptr(c.mloc()*n),&wftmp[0]);
      // pointers to rmat
      vector<double *> p_rmat( ndiv_[0]*ndiv_[1] );
      for (int iproj=0; iproj<ndiv_[0]*ndiv_[1]; iproj++)
      {
        int index = n * rmat_[iproj]->mloc();
        p_rmat[iproj] = rmat_[iproj]->valptr(index);
      }
      // copy wf to rmat arrays
      for ( int iz = 0; iz < np2loc_; iz++ )
      {
        for ( int ixy = 0; ixy < np01_; ixy++ )
        {
          int xy_proj = xy_proj_[ixy];
          (*p_rmat[xy_proj]) = wftmp[ixy + iz*np01_].real();
          p_rmat[xy_proj]++;
        }
      }
    }
#ifdef TIMING
    tmap["wf FFT"].stop();
#endif

    // compute the x/y A matrices
#ifdef TIMING
      tmap["xy products"].start();
#endif
    {
      // clear a matrices
      for ( int i = 0; i < nmat_; i++ )
      {
        amat_[i]->clear();
      }

      int size=amat_[0]->size();

      // allocate matrix for products of projected parts
      DoubleMatrix products(c.context(),c.n(),c.n(),c.nb(),c.nb());

      for ( int i_proj=0; i_proj<rmat_.size(); i_proj++ )
      {
        // get index of the subbox
        int i_slice[2];
        i_slice[1]= i_proj/ndiv_[0];
        i_slice[0]= i_proj-ndiv_[0]*i_slice[1];

        // compute product of projections in real space
        products.gemm('t','n',1.0,(*rmat_[i_proj]),(*rmat_[i_proj]),0.0);

        // add product to the x/y A matrices
        for ( int ilevel = 0, imat = 0; ilevel < nlevelsmax_; ilevel++ )
        {
          for ( int idir = 0; idir < 3; idir++ )
          {
            // matrix A is an x or y matrix
            if ( ilevel<nlevels_[idir] && idir<2 )
            {
              if ( walsh( ilevel , ndiv_[idir] , i_slice[idir] ) )
              {
                // add product to the matrix
                double *coeff_source=products.valptr(0);
                double *coeff_destination=amat_[imat]->valptr(0);
                for ( int i=0; i<size; i++ )
                  coeff_destination[i]+=coeff_source[i];
              }
              imat++;
            }
            // else: matrix A is a z matrix
            else if ( ilevel<nlevels_[idir] )
            {
              imat++;
            }
          }
        }
      }

      // normalize xy matrices
      for ( int ilevel = 0, imat = 0; ilevel<nlevelsmax_; ilevel++ )
      {
        for ( int idir = 0; idir < 3; idir++ )
        {
          // matrix A is a projector in the  x or y direction
          if ( ilevel < nlevels_[idir] && idir<2 )
          {
            // normalize coeffs
            double *coeff=amat_[imat]->valptr(0);
            double fac = 1.0 / (np_[0]*np_[1]*np_[2]);
            for ( int i = 0; i < size; i++ )
              coeff[i] *= fac;
            imat++;
          }
          // else: matrix A is a projector in the z direction
          else if ( ilevel<nlevels_[idir] )
          {
            imat++;
          }
        }
      }
    }
#ifdef TIMING
    tmap["xy products"].stop();
#endif

    // compute the z A matrices
#ifdef TIMING
    tmap["z products"].start();
#endif
    {
      int imat=0;
      // matrix for projected wf
      ComplexMatrix c_pz(c.context(),c.m(),c.n(),c.mb(),c.nb());;
      // proxy for complex->real matrix product
      DoubleMatrix c_proxy(c);
      DoubleMatrix c_pz_proxy(c_pz);

      for ( int ilevel=0; ilevel<nlevelsmax_; ilevel++ )
      {
        // loop on directions
        for ( int idir=0; idir<3; idir++ )
        {
          // matrix A is for x or y direction
          if ( ilevel<nlevels_[idir] && idir<2 )
          {
            imat++;
          }
          // else: matrix A is for z direction
          else if ( ilevel<nlevels_[idir] )
          {
            for ( int n = 0; n < nstloc; n++ )
            {
              // p_rmat: pointers to rmat arrays
              vector<double *> p_rmat( ndiv_[0]*ndiv_[1] );
              for ( int iproj=0; iproj<ndiv_[0]*ndiv_[1]; iproj++ )
              {
                int index = n * rmat_[iproj]->mloc();
                p_rmat[iproj] = rmat_[iproj]->valptr(index);
              }
              // save values of wf*walsh in rmat arrays
              for ( int iz = 0; iz < np2loc_; iz++ )
              {
                int izglobal = ft_->np2_first() + iz;
                double walsh_z = walsh(ilevel, np_[2], izglobal);
                for ( int ixy = 0; ixy < np01_; ixy++ )
                {
                  int i = ixy + iz * np01_;
                  int xy_proj = xy_proj_[ixy];
                  wftmp[i] = (*p_rmat[xy_proj]) * walsh_z;
                  p_rmat[xy_proj]++;
                }
              }
              ft_->forward(&wftmp[0],c_pz.valptr(c_pz.mloc()*n));
            }
            // compute the product
            // factor -2.0 in next line: G and -G
            amat_[imat]->gemm('t','n',2.0,c_pz_proxy,c_proxy,0.0);
            // rank-1 update using first row of cd_proxy() and c_proxy
            amat_[imat]->ger(-1.0,c_pz_proxy,0,c_proxy,0);
            imat++;
          }
        }
      }
    }
#ifdef TIMING
    tmap["z products"].stop();
#endif

#ifdef DEBUG
    // check the values of the amat matrices
#ifdef TIMING
    tmap["check_amat"].start();
#endif
    check_amat(c);
#ifdef TIMING
    tmap["check_amat"].stop();
#endif
#endif

    // set to zero matrix elements of the matrices amat_[i] if they couple
    // states with differing occupation numbers
    trim_amat(occ);

#ifdef DEBUG_PRINT_MAT
    for ( int k = 0; k < amat_.size(); k++ )
    {
      cout << "A(k=" << k << "):" << endl;
      cout << *amat_[k];
    }
#endif

    // joint approximate diagonalization of the matrices amat (jade)
#ifdef TIMING
    tmap["jade"].start();
#endif

    // diagonal values adiag[k][i]
    vector<vector<double> > adiag(nmat_);

    // diagonalize projectors
    int maxsweep=200;
    double tol=max(1E-10,epsilon/10);
    jade(maxsweep,tol,amat_,*u_,adiag);

#ifdef TIMING
    tmap["jade"].stop();
#endif

#ifdef DEBUG_PRINT_MAT
    cout << "U:" << endl;
    cout << *u_;
#endif

    // compute localization vector
#ifdef TIMING
    tmap["localization"].start();
#endif

    for ( int n = 0; n < nst; n++ )
      localization_[ispin][n] = 0;

    for ( int imat=0; imat<nmat_; imat++ )
    {
      // broadcast diagonal of all matrices a to all tasks
      MPI_Bcast( (void *) &adiag[imat][0], adiag[imat].size(),
        MPI_DOUBLE, 0, wf.context().comm() );

      for ( int n = 0; n < nst; n++ )
      {
        if ( adiag[imat][n] < epsilon )
        {
          localization_[ispin][n] += 1<<(2*imat);
        }
        else if ( adiag[imat][n] > 1.0-epsilon )
        {
          localization_[ispin][n] += 1<<(2*imat+1);
        }
        else
        {
          localization_[ispin][n] += (1<<(2*imat)) + (1<<(2*imat+1));
        }
      }
    }
#ifdef TIMING
    tmap["localization"].stop();
#endif

#ifdef DEBUG
    // print localization vector and number of overlaps (including self)
    // for each state
    if ( gcontext_.onpe0() )
    {
      int sum = 0;
      for ( int i = 0; i < nst; i++ )
      {
        int count = 0;
        for ( int j = 0; j < nst; j++ )
        {
          if ( overlap(ispin,i,j) )
            count++;
        }
        cout << "localization[" << i << "]: "
             << localization_[ispin][i] << " "
             << bitset<30>(localization_[ispin][i]) << "  overlaps: "
             << count << endl;
        sum += count;
      }
      cout << "total overlaps: " << sum << " / " << nst*nst
           << " = " << ((double) sum)/(nst*nst) << endl;
    }
#endif

    // broadcast localization to all tasks to ensure consistency
    MPI_Bcast( (void *) &localization_[ispin][0], localization_[ispin].size(),
      MPI_LONG, 0, wf.context().comm() );

    // distribute the states among the process columns
#ifdef TIMING
    tmap["distribute"].start();
#endif

    // distribute only if the keyword BISECTION_NODIST is not in the debug var
    if ( s_.ctrl.debug.find("BISECTION_NODIST") == string::npos )
    {
      distribute(ispin);
    }
    else
    {
      if ( gcontext_.onpe0() )
        cout << " Bisection::localize: distribution disabled" << endl;
    }

    // distribute only if the keyword BISECTION_NODIST is not in the debug var
    if ( s_.ctrl.debug.find("BISECTION_NODIST") == string::npos )
    {
      distribute(ispin);
    }
    else
    {
      if ( gcontext_.onpe0() )
        cout << " Bisection::localize: distribution disabled" << endl;
    }

#if DEBUG
    // print localization vectors and overlaps after distribution
    if ( gcontext_.onpe0() )
    {
      int sum = 0;
      for ( int i = 0; i < nst; i++ )
      {
        int count = 0;
        for ( int j = 0; j < nst; j++ )
        {
          if ( overlap(ispin,i,j) )
            count++;
        }
        cout << "localization[" << i << "]: "
             << localization_[ispin][i] << " "
             << bitset<30>(localization_[ispin][i]) << "  overlaps: "
             << count << endl;
        sum += count;
      }
      cout << " Bisection::localize: total overlaps: " << sum << " / "
           << nst*nst << " = " << ((double) sum)/(nst*nst) << endl;
    }
#endif

#ifdef TIMING
    tmap["distribute"].stop();
#endif

    // apply rotation to the states
#ifdef TIMING
    tmap["rotate"].start();
#endif
    ComplexMatrix cp(c);
    DoubleMatrix cp_proxy(cp);
    DoubleMatrix c_proxy(c);
    c_proxy.gemm('n','n',1.0,cp_proxy,*u_,0.0);
#ifdef TIMING
    tmap["rotate"].stop();
#endif

  } // ispin

  // print timers
#ifdef TIMING
  for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    wf.context().dmin(1,1,&tmin,1);
    wf.context().dmax(1,1,&tmax,1);
    if ( wf.context().myproc()==0 )
    {
      cout << "<timing name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
    }
  }
#endif
}

#define DEGREE_ORDERING
////////////////////////////////////////////////////////////////////////////////
#ifndef  DEGREE_ORDERING
struct LocalizationLess
{
  // Function object to compare localization vectors
  public:
  vector<long int>& loc_;
  LocalizationLess(vector<long int>& loc) : loc_(loc) {}
  bool operator() (int i, int j) const
//  { return comesBefore(loc_[i],loc_[j]); }
// bool comesBefore(long int loc1, long int loc2)
{
  long int loc1 = loc_[i];
  long int loc2 = loc_[j];

  // fast return if possible
  if ( loc1==loc2 ) return true;

  // count the number of non zero projections
  int n_proj1=0;
  int n_proj2=0;
  {
    long int loc1_=loc1;
    long int loc2_=loc2;

    while ( loc1_!=0 && loc2_!=0 )
    {
      // right projection
      n_proj1 += loc1_ & 1;
      n_proj2 += loc2_ & 1;
      loc1_ >>= 1;
      loc2_ >>= 1;
      // left projection
      n_proj1 += loc1_ & 1;
      n_proj2 += loc2_ & 1;
      loc1_ >>= 1;
      loc2_ >>= 1;
    }
  }

  // if the number of projections are different
  if      ( n_proj1>n_proj2 ) return true;
  else if ( n_proj1<n_proj2 ) return false;

  // if the number of projections are the same,
  // compare on which is localized first
  if (0)
  {
    long int loc1_=loc1;
    long int loc2_=loc2;

    while ( loc1_!=0 && loc2_!=0 )
    {
      // get the weight of projections for the two states
      bool p_right1 = loc1_ & 1;
      bool p_left1  = loc1_ & 2;
      bool p_right2 = loc2_ & 1;
      bool p_left2  = loc2_ & 2;

      // compare projections
      if ( p_right1 < p_right2 || p_left1 > p_left2 ) return false;
      if ( p_right1 > p_right2 || p_left1 < p_left2 ) return true;

      // shift bitset
      loc1_ >>= 2;
      loc2_ >>= 2;
    }
  }
  else
  {
    long int loc1_=loc1;
    long int loc2_=loc2;

    while ( loc1_!=0 && loc2_!=0 )
    {
      // get the weight of projections for the two states
      // on the right projectors
      bool p_right1 = loc1_ & 1;
      bool p_right2 = loc2_ & 1;

      // compare projections
      if ( p_right1 < p_right2 ) return false;
      if ( p_right1 > p_right2 ) return true;

      // shift bitset
      loc1_ >>= 2;
      loc2_ >>= 2;
    }
    // reset localizations and shift for left projection analysis
    loc1_= loc1 >> 1;
    loc2_= loc2 >> 1;
    //
    while ( loc1_!=0 && loc2_!=0 )
    {
      // get the weight of projections for the two states
      // on the right projectors
      bool p_left1 = loc1_ & 1;
      bool p_left2 = loc2_ & 1;

      // compare projections
      if ( p_left1 > p_left2 ) return false;
      if ( p_left1 < p_left2 ) return true;

      // shift bitset
      loc1_ >>= 2;
      loc2_ >>= 2;
    }
  }

  // should never get here
  cout << "Error: should not reach this point in function "
       << "- comesBefore(long int loc1, long int loc2) -" << endl;
  cout << "loc1=" << loc1 << endl;
  cout << "loc2=" << loc2 << endl;

  exit(0);
}
}; // end of struct LocalizationLess
#endif

////////////////////////////////////////////////////////////////////////////////
void Bisection::distribute(int ispin)
{
  // permute the columns of u according to the order defined by the
  // localization vector
  vector<long int>& localization = localization_[ispin];

  vector<int> index(localization.size());
  for ( int j = 0; j < index.size(); j++ )
    index[j] = j;

  // compute the degree of the vertices of the exchange graph
  vector<int> degree(localization.size());
  for ( int i = 0; i < localization.size(); i++ )
  {
    int count = 0;
    for ( int j = 0; j < localization.size(); j++ )
    {
      if ( overlap(ispin,i,j) )
        count++;
    }
    degree[i] = count;
  }

#ifdef DEGREE_ORDERING
  // Create function object for comparison of degree
  VectorLess<int> degree_less(degree);
  sort(index.begin(), index.end(), degree_less);
  // At this point degree[index[i]] <= degree[index[j]] if i < j
  for ( int i = 0; i < index.size()-1; i++ )
    assert(degree[index[i]] <= degree[index[i+1]]);
#else
  // Create function object for comparison of localization vectors
  LocalizationLess loc_less(localization);
  sort(index.begin(), index.end(), loc_less);
#endif

#if DEBUG
  if ( u_->context().onpe0() )
  {
    cout << "degree order after sort:" << endl;
    for ( int j = 0; j < index.size(); j++ )
      cout << j << " -> " << index[j]
           << "  " << degree[index[j]]
           << endl;
  }
#endif

  // distribute the states to process columns in round robin fashion
  // Assume that the states are initially ordered by increasing degree
  // i.e. degree(index[i]) < degree(index[j]) if i < j

  const int n = index.size();
  const int nb = u_->nb();

  vector<int> distrib_index(n);
  int ibase = 0;
  int icol = 0;
  for ( int i = 0; i < n; i++ )
  {
    // check if next slot is beyond n
    if ( ibase + icol * nb >= n )
    {
      // restart in column 0 with incremented ibase
      icol = 0;
      ibase++;
    }
    distrib_index[ibase + icol * nb] = i;
    icol++;
  }

  // combine index[i] and distrib_index[i]
  vector<int> itmp(index.size());
  for ( int i = 0; i < n; i++ )
    itmp[i] = index[distrib_index[i]];
  index = itmp;

  // apply the permutation defined by index to localization
  vector<long int> loc_tmp(localization.size());
  for ( int i = 0; i < index.size(); i++ )
  {
    loc_tmp[i] = localization[index[i]];
  }
  localization = loc_tmp;

  // apply the permutation defined by index to occupation numbers
  vector<double> occ_tmp(s_.wf.sd(ispin,0)->nst());
  const double* occ = s_.wf.sd(ispin,0)->occ_ptr();
  for ( int i = 0; i < index.size(); i++ )
  {
    occ_tmp[i] = occ[index[i]];
  }
  s_.wf.sd(ispin,0)->set_occ(occ_tmp);

  // transform the permutation defined by index[i] into a sequence of
  // transpositions stored in the vector pivot[i]
  vector<int> pivot;
  for ( int i=0; i<index.size(); i++ )
  {
    for ( int j=i; j<index.size(); j++ )
    {
      if ( index[j] == i )
      {
        int tmp = index[i];
        index[i] = index[j];
        index[j] = tmp;
        pivot.push_back(j);
        break;
      }
    }
  }

  // create a local pivot vector on this process (size u.nloc())
  // this vector must be replicated on all tasks of the process grid columns
  vector<int> locpivot(u_->nloc());
  // fill the local pivot vector on all tasks
  // add 1 to index values for lapack fortran index convention (start at 1)
  for ( int j=0; j < u_->nloc(); j++ )
  {
    int jglobal = u_->jglobal(j);
    locpivot[j] = pivot[jglobal]+1;
  }

  // apply the permutation
  u_->lapiv('B','C',&locpivot[0]);

#if DEBUG
  // recompute the degree of the vertices of the exchange graph
  for ( int i = 0; i < localization.size(); i++ )
  {
    int count = 0;
    for ( int j = 0; j < localization.size(); j++ )
    {
      if ( overlap(ispin,i,j) )
        count++;
    }
    degree[i] = count;
  }

  if ( u_->context().onpe0() )
  {
    cout << "degree after permutation:" << endl;
    for ( int j = 0; j < degree.size(); j++ )
      cout << " deg[" << j << "] = " << degree[j] << endl;

    cout << "occ after permutation:" << endl;
    for ( int j = 0; j < degree.size(); j++ )
      cout << " occ[" << j << "] = " << s_.wf.sd(ispin,0)->occ(j) << endl;
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
bool Bisection::check_amat(ComplexMatrix &c)
{
  // allocate memory for copy of the wave functions
  ComplexMatrix cd(c.context(),c.m(),c.n(),c.mb(),c.nb());

  // precompute values of Walsh functions in three directions
  // wx[l][i] = value of level l Walsh function at position i in direction x
  vector<vector<vector<int> > > w(3);
  for ( int idir=0; idir<3; idir++ )
  {
    w[idir].resize(nlevels_[idir]);
    //
    for ( int l = 0; l < nlevels_[idir]; l++ )
    {
      w[idir][l].resize(np_[idir]);
      for ( int i = 0; i < np_[idir]; i++ )
      {
        w[idir][l][i] = walsh(l,np_[idir],i);
      }
    }
  }

  // compute matrices B_k = <wf|dwf>
  vector<DoubleMatrix*> bmat(nmat_);
  for ( int k = 0; k < bmat.size(); k++ )
  {
    bmat[k] = new DoubleMatrix(c.context(),c.n(),c.n(),c.nb(),c.nb());
  }
  DoubleMatrix c_proxy(c);
  DoubleMatrix cd_proxy(cd);

  vector<complex<double> > wftmp(ft_->np012loc());
  complex<double> *f = &wftmp[0];

  // compute matrices A at all levels
  for ( int l=0 , imat=0; l < nlevelsmax_; l++ )
  {
    // x direction
    if ( l<nlevels_[0] )
    {
      cd_proxy = c_proxy;
      for ( int n = 0; n < c.nloc(); n++ )
      {
        for ( int i = 0; i < np012loc_; i++ )
          f[i] = 0.0;
        ft_->backward(cd.cvalptr(cd.mloc()*n),&wftmp[0]);
        for ( int i = 0; i < np_[0]; i++ )
        {
          if ( w[0][l][i] == 0 )
          {
            for ( int j = 0; j < np_[1]; j++ )
              for ( int k = 0; k < ft_->np2_loc(); k++ )
                f[i +  np_[0] * ( j + np_[1] * k )] = 0.0;
          }
        }
        ft_->forward(&wftmp[0],cd.valptr(cd.mloc()*n));
      }
      // factor -2.0 in next line: G and -G
      bmat[imat]->gemm('t','n',2.0,cd_proxy,c_proxy,0.0);
      // rank-1 update using first row of cd_proxy() and c_proxy
      bmat[imat]->ger(-1.0,cd_proxy,0,c_proxy,0);
      imat++;
    }

    // y direction
    if ( l<nlevels_[1] )
    {
      cd_proxy = c_proxy;
      for ( int n = 0; n < c.nloc(); n++ )
      {
        for ( int i = 0; i < np012loc_; i++ )
          f[i] = 0.0;
        ft_->backward(cd.cvalptr(cd.mloc()*n),&wftmp[0]);
        for ( int j = 0; j < np_[1]; j++ )
        {
          if ( w[1][l][j] == 0 )
          {
            for ( int i = 0; i < np_[0]; i++ )
              for ( int k = 0; k < ft_->np2_loc(); k++ )
                f[i +  np_[0] * ( j + np_[1] * k )] = 0.0;
          }
        }
        ft_->forward(&wftmp[0],cd.valptr(cd.mloc()*n));
      }
      // factor -2.0 in next line: G and -G
      bmat[imat]->gemm('t','n',2.0,cd_proxy,c_proxy,0.0);
      // rank-1 update using first row of cd_proxy() and c_proxy
      bmat[imat]->ger(-1.0,cd_proxy,0,c_proxy,0);
      imat++;
    }

    // z direction
    if ( l<nlevels_[2] )
    {
      cd_proxy = c_proxy;
      for ( int n = 0; n < c.nloc(); n++ )
      {
        for ( int i = 0; i < np012loc_; i++ )
          f[i] = 0.0;
        ft_->backward(cd.cvalptr(cd.mloc()*n),&wftmp[0]);
        for ( int k = 0; k < ft_->np2_loc(); k++ )
        {
          int kglobal = ft_->np2_first() + k;
          const int istart = k * np_[0] * np_[1];
          if ( w[2][l][kglobal] == 0 )
          {
            for ( int ij = 0; ij < np_[0]*np_[1]; ij++ )
              f[istart+ij] = 0.0;
          }
        }
        ft_->forward(&wftmp[0],cd.valptr(cd.mloc()*n));
      }
      // factor -2.0 in next line: G and -G
      bmat[imat]->gemm('t','n',2.0,cd_proxy,c_proxy,0.0);
      // rank-1 update using first row of cd_proxy() and c_proxy
      bmat[imat]->ger(-1.0,cd_proxy,0,c_proxy,0);
      imat++;
    }
  } // for l

  // testing the matrices
  for ( int imat=0; imat<nmat_; imat++ )
  {
    double *a=amat_[imat]->valptr(0);
    double *b=bmat[imat]->valptr(0);
    int ncoeff=amat_[imat]->mloc()*amat_[imat]->nloc();

    for ( int i=0; i<ncoeff; i++ )
    {
      if (fabs(a[i]-b[i])>1e-10)
      {
        cout << "error > 1.e-10 for matrix " << imat << endl;
      }
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::trim_amat(const vector<double>& occ)
{
  // set to zero the matrix elements of the matrices amat_[k] if they couple
  // states with differing occupation numbers

  // check if all occupation numbers are the same
  double occ_max = 0.0, occ_min = 2.0;
  for ( int i = 0; i < occ.size(); i++ )
  {
    occ_max = max(occ_max, occ[i]);
    occ_min = min(occ_min, occ[i]);
  }
  // return if all occupation numbers are equal
  if ( fabs(occ_max-occ_min) < 1.e-12 )
    return;

  const int mloc = amat_[0]->mloc();
  const int nloc = amat_[0]->nloc();
  // loop over elements local to this task
  for ( int i = 0; i < mloc; i++ )
  {
    const int iglobal = amat_[0]->iglobal(i);
    for ( int j = 0; j < nloc; j++ )
    {
      const int jglobal = amat_[0]->jglobal(j);

      const int ival = i + mloc * j;
      if ( fabs(occ[iglobal] - occ[jglobal]) > 1.e-12 )
        for ( int k = 0; k < amat_.size(); k++ )
          (*amat_[k])[ival] = 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// overlap: return true if the functions i and j overlap
bool Bisection::overlap(int ispin, int i, int j)
{
  long int localization_i = localization_[ispin][i];
  long int localization_j = localization_[ispin][j];
  while ( localization_i!=0 && localization_j!=0 )
  {
    // get the weight of projections for each state
    bool p_right_i = localization_i & 1;
    bool p_left_i  = localization_i & 2;
    bool p_right_j = localization_j & 1;
    bool p_left_j  = localization_j & 2;

    // return false as soon as the states are found to be separated
    if ( !( ( p_right_i && p_right_j ) || ( p_left_i && p_left_j ) ) )
      return false;

    localization_i >>= 2;
    localization_j >>= 2;
  }

  // return true if the states overlap
  return true;
}

////////////////////////////////////////////////////////////////////////////////
double Bisection::pair_fraction(int ispin)
{
  // pair_fraction: return fraction of pairs having non-zero overlap
  const int nst = s_.wf.nst(ispin);
  int sum = 0;
  for ( int i = 0; i < nst; i++ )
  {
    int count = 0;
    for ( int j = 0; j < nst; j++ )
    {
      if ( overlap(ispin,i,j) )
        count++;
    }
    sum += count;
  }
  return ((double) sum)/(nst*nst);
}

////////////////////////////////////////////////////////////////////////////////
double Bisection::size(int ispin, int i)
{
  // size: return fraction of the domain on which state i is non-zero
  long int loc = localization_[ispin][i];
  double size = 1.0;
  // process all projectors
  while ( loc != 0 )
  {
    // weight of projections
    bool p_right = loc & 1;
    bool p_left  = loc & 2;

    if ( (p_right && !p_left) || (!p_right && p_left) )
      size *= 0.5;

    loc >>= 2;
  }

  // return true if the states overlap
  return size;
}

////////////////////////////////////////////////////////////////////////////////
double Bisection::total_size(int ispin)
{
  // total_size: return normalized sum of sizes
  const int nst = s_.wf.nst(ispin);
  double sum = 0.0;
  for ( int i = 0; i < nst; i++ )
    sum += size(ispin,i);
  return sum / nst;
}
