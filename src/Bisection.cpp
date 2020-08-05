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
// Bisection.cpp
//
////////////////////////////////////////////////////////////////////////////////
#include "Bisection.h"
#include <bitset>
#include <algorithm>
#include "jade.h"
#include "FourierTransform.h"

using namespace std;

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
Bisection::Bisection(const SlaterDet& sd, const int nlevels[3])
  : ctxt_(sd.context())
{
  // localization vectors are long int
  assert(sizeof(long int) >= 4);

  nst_ = sd.nst();
  nstloc_ = sd.nstloc();

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

  assert( nlevelsmax_ >0 );

  // real-space grid size for wave functions
  np_[0] = sd.basis().np(0);
  np_[1] = sd.basis().np(1);
  np_[2] = sd.basis().np(2);
  // adapt the grid dimensions to the levels of bisection
  for ( int idim = 0; idim < 3; idim++ )
  {
    for ( int j=1; j<=nlevels[idim]; j++ )
    {
      int base = 1 << j;
      if ( np_[idim] % base != 0 ) np_[idim] += base/2;
    }
  }
  while (!sd.basis().factorizable(np_[0])) np_[0] += (1<<nlevels[0]);
  while (!sd.basis().factorizable(np_[1])) np_[1] += (1<<nlevels[1]);
  while (!sd.basis().factorizable(np_[2])) np_[2] += (1<<nlevels[2]);

  // number of grid points of augmented grid for normalization
  ft_ = new FourierTransform(sd.basis(),np_[0],np_[1],np_[2]);
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

  // each projector uses two bits of the localization vector
  // check that sizeof(long int) is sufficient
  // number of bits in a long int is 8*sizeof(long int)
  assert(2*nmat_ <= 8*sizeof(long int));

  // allocate A matrices
  amat_.resize(nmat_);
  const ComplexMatrix &c = sd.c();
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
      (int)xyproj_rsize.size(), MPI_INT ,MPI_MAX , ctxt_.comm());

    // allocate matrices rmat_[i]
    for ( int i = 0; i < rmat_.size(); i++ )
    {
      int n = c.n();
      int nb = c.nb();
      int m = max_xyproj_rsize[i] * c.context().nprow();
      int mb = max_xyproj_rsize[i];
      rmat_[i] = new DoubleMatrix(c.context(),m,n,mb,nb);
      rmat_[i]->clear();
    }
  }

  localization_.resize(nst_);
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
void Bisection::compute_transform(const SlaterDet& sd, int maxsweep, double tol)
{
  // compute the transformation with tolerance tol

  map<string,Timer> tmap;

  // check that basis is real
  // jade is not implemented for complex matrices
  assert( sd.basis().real() );

  const ComplexMatrix& c = sd.c();
  const vector<double>& occ = sd.occ();
  assert(occ.size() == nst_);

#ifdef TIMING
  tmap["wf FFT"].start();
#endif
  // Fourier transform states and save real-space values in
  // rmat_ matrices
  vector<complex<double> > wftmp(np012loc_,0.0);
  for ( int n = 0; n < nstloc_; n++ )
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
  // clear a matrices
  for ( int i = 0; i < nmat_; i++ )
    amat_[i]->clear();

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
            const int size=amat_[imat]->size();
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
        const int size=amat_[imat]->size();
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

#ifdef TIMING
  tmap["xy products"].stop();
#endif

  // compute the z A matrices
#ifdef TIMING
  tmap["z products"].start();
#endif
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
        for ( int n = 0; n < nstloc_; n++ )
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

  // diagonal values adiag_[k][i]
  // adiag_ is resized by jade

  // diagonalize projectors
  // int nsweep = jade(maxsweep,tol,amat_,*u_,adiag_);
  jade(maxsweep,tol,amat_,*u_,adiag_);
#ifdef TIMING
  if ( ctxt_.onpe0() )
    cout << "Bisection::compute_transform:"
         << " maxsweep=" << maxsweep << " tol=" << tol << endl;
#endif

#ifdef TIMING
  tmap["jade"].stop();
#endif

#ifdef DEBUG_PRINT_MAT
  cout << "U:" << endl;
  cout << *u_;
#endif

  for ( int imat=0; imat<nmat_; imat++ )
  {
    // broadcast diagonal of all matrices a to all tasks
    MPI_Bcast( (void *) &adiag_[imat][0], adiag_[imat].size(),
               MPI_DOUBLE, 0, ctxt_.comm() );
  }
  // print timers
#ifdef TIMING
  for ( map<string,Timer>::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    double sbuf = tmin;
    double rbuf = 0.0;
    MPI_Reduce(&sbuf,&rbuf,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    sbuf = tmax;
    rbuf = 0.0;
    MPI_Reduce(&sbuf,&rbuf,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( MPIdata::onpe0() )
    {
      string s = "name=\"" + (*i).first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::compute_localization(double epsilon)
{
  // compute localization vector for a threshold epsilon
  for ( int n = 0; n < nst_; n++ )
  {
    localization_[n] = 0;
    for ( int imat = 0; imat < nmat_; imat++ )
    {
      if ( adiag_[imat][n] < epsilon )
        localization_[n] += 1<<(2*imat);
      else if ( adiag_[imat][n] > 1.0-epsilon)
        localization_[n] += 1<<(2*imat+1);
      else
        localization_[n] += (1<<(2*imat)) + (1<<(2*imat+1));
    }
  }

#ifdef DEBUG
  // print localization vector and number of overlaps (including self)
  // for each state
  if ( ctxt_.onpe0() )
  {
    int sum = 0;
    for ( int i = 0; i < nst_; i++ )
    {
      int count = 0;
      for ( int j = 0; j < nst_; j++ )
      {
        if ( overlap(i,j) )
          count++;
      }
      cout << "localization[" << i << "]: "
           << localization_[i] << " "
           << bitset<30>(localization_[i]) << "  overlaps: "
           << count << endl;
      sum += count;
    }
    cout << "total overlaps: " << sum << " / " << nst_*nst_
         << " = " << ((double) sum)/(nst_*nst_) << endl;
  }
#endif

  // broadcast localization to all tasks to ensure consistency
  MPI_Bcast( (void *) &localization_[0], localization_.size(),
             MPI_LONG, 0, ctxt_.comm() );

}

////////////////////////////////////////////////////////////////////////////////
void Bisection::forward(SlaterDet& sd)
{
  forward(*u_,sd);
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::forward(DoubleMatrix& u, SlaterDet& sd)
{
  // apply the bisection transformation to the SlaterDet sd
  // apply the rotation u to sd.c()
  ComplexMatrix& c = sd.c();
  ComplexMatrix cp(c);
  DoubleMatrix cp_proxy(cp);
  DoubleMatrix c_proxy(c);
  c_proxy.gemm('n','n',1.0,cp_proxy,u,0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::backward(SlaterDet& sd)
{
  backward(*u_,sd);
}

////////////////////////////////////////////////////////////////////////////////
void Bisection::backward(DoubleMatrix& u, SlaterDet& sd)
{
  // apply the inverse bisection transformation to SlaterDet sd
  // apply rotation u^T to sd
  ComplexMatrix& c = sd.c();
  ComplexMatrix cp(c);
  DoubleMatrix cp_proxy(cp);
  DoubleMatrix c_proxy(c);
  c_proxy.gemm('n','t',1.0,cp_proxy,u,0.0);
}

////////////////////////////////////////////////////////////////////////////////
bool Bisection::check_amat(const ComplexMatrix &c)
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

  const double trim_tol = 1.e-6;
  // check if all occupation numbers are the same
  double occ_max = 0.0, occ_min = 2.0;
  for ( int i = 0; i < occ.size(); i++ )
  {
    occ_max = max(occ_max, occ[i]);
    occ_min = min(occ_min, occ[i]);
  }
  // return if all occupation numbers are equal
  if ( fabs(occ_max-occ_min) < trim_tol )
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
      if ( fabs(occ[iglobal] - occ[jglobal]) > trim_tol )
        for ( int k = 0; k < amat_.size(); k++ )
          (*amat_[k])[ival] = 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
bool Bisection::overlap(int i, int j) const
{
  return overlap(localization_,i,j);
}

////////////////////////////////////////////////////////////////////////////////
bool Bisection::overlap(const vector<long int>& loc_, int i, int j) const
{
  // overlap: return true if the functions i and j overlap according
  // to the localization vector loc_
  long int loc_i = loc_[i];
  long int loc_j = loc_[j];
  while ( loc_i!=0 && loc_j!=0 )
  {
    // get the weight of projections for each state
    bool p_right_i = loc_i & 1;
    bool p_left_i  = loc_i & 2;
    bool p_right_j = loc_j & 1;
    bool p_left_j  = loc_j & 2;

    // return false as soon as the states are found to be separated
    if ( !( ( p_right_i && p_right_j ) || ( p_left_i && p_left_j ) ) )
      return false;

    loc_i >>= 2;
    loc_j >>= 2;
  }

  // return true if the states overlap
  return true;
}

////////////////////////////////////////////////////////////////////////////////
double Bisection::pair_fraction(void) const
{
  // pair_fraction: return fraction of pairs having non-zero overlap
  // count pairs (i,j) having non-zero overlap for i != j only
  int sum = 0;
  for ( int i = 1; i < nst_; i++ )
  {
    int count = 0;
    for ( int j = i+1; j < nst_; j++ )
    {
      if ( overlap(i,j) )
        count++;
    }
    sum += count;
  }
  // add overlap with self: (i,i)
  sum += nst_;
  return ((double) sum)/((nst_*(nst_+1))/2);
}

////////////////////////////////////////////////////////////////////////////////
double Bisection::size(int i) const
{
  // size: return fraction of the domain on which state i is non-zero
  long int loc = localization_[i];
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
double Bisection::total_size(void) const
{
  // total_size: return normalized sum of sizes
  double sum = 0.0;
  for ( int i = 0; i < nst_; i++ )
    sum += size(i);
  return sum / nst_;
}
