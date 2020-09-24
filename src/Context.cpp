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
// Context.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Context.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include "blacs.h"
#include <mpi.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Context::Context(MPI_Comm comm) : ictxt_(-1), myrow_(-1), mycol_(-1)
{
  // construct a single-row Context
  int nprocs;
  char order='R';
  comm_ = comm;
  MPI_Comm_size(comm_,&nprocs);
  MPI_Comm_rank(comm_,&mype_);
  ictxt_ = Csys2blacs_handle(comm_);
  nprow_ = 1;
  npcol_ = nprocs;

  Cblacs_gridinit( &ictxt_, &order, nprow_, npcol_ );

  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);

  size_ = nprow_ * npcol_;
  myproc_ = myrow_ < 0 ? -1 : mycol_ + npcol_ * myrow_;
  onpe0_ = ( mype_ == 0 );
  active_ = ( ictxt_ >= 0 );

  pmap_.resize(size_);
  for ( int i = 0; i < size_; i++ )
    pmap_[i] = i;
}

////////////////////////////////////////////////////////////////////////////////
Context::Context(MPI_Comm comm, int nprow, int npcol) :
  ictxt_(-1), myrow_(-1), mycol_(-1), nprow_(nprow), npcol_(npcol)
{
  int nprocs;
  char order = 'C';
  comm_ = comm;
  MPI_Comm_size(comm_,&nprocs);
  MPI_Comm_rank(comm_,&mype_);
  ictxt_ = Csys2blacs_handle(comm_);
  if ( nprocs < nprow * npcol )
  {
    cout << " nprocs=" << nprocs << endl;
    cout << " Context nprow*npcol > nprocs" << endl;
    Cblacs_abort(ictxt_, 1);
  }
  Cblacs_gridinit( &ictxt_, &order, nprow, npcol );

  // get values of nprow_, npcol_, myrow_ and mycol_ in the new context
  if ( ictxt_ >= 0 )
    Cblacs_gridinfo(ictxt_, &nprow_, &npcol_, &myrow_, &mycol_);

  size_ = nprow_ * npcol_;
  myproc_ = Cblacs_pnum(ictxt_,myrow_,mycol_);
  onpe0_ = ( mype_ == 0 );
  active_ = ( ictxt_ >= 0 );

  pmap_.resize(size_);
  // column-major order
  int i = 0;
  for ( int ic = 0; ic < npcol; ic++ )
    for ( int ir = 0; ir < nprow; ir++ )
    {
      pmap_[ir+nprow*ic] = i;
      i++;
    }
}

////////////////////////////////////////////////////////////////////////////////
Context::~Context()
{
  if ( myrow_ != -1 )
  {
    // cout << " Context destructor called on ictxt = " << ictxt_ << endl;
    Cblacs_gridexit( ictxt_ );
  }
}

////////////////////////////////////////////////////////////////////////////////
void Context::dsend(int m, int n, double* a, int lda, int rdest, int cdest)
  const { Cdgesd2d(ictxt_,m,n,a,lda,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::drecv(int m, int n, double* a, int lda, int rsrc, int csrc)
  const { Cdgerv2d(ictxt_,m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::dsum(char scope, char topology, int m, int n, double* a,
  int lda, int rdest, int cdest) const
{ Cdgsum2d(ictxt_,&scope,&topology,m,n,a,lda,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::dsum(char scope, int m, int n, double* a, int lda) const
{ dsum(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dsum(int m, int n, double* a, int lda) const
{ dsum('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmax(char scope, char topology, int m, int n, double* a,
  int lda, int* ra, int* ca, int rcflag, int rdest, int cdest) const
{ Cdgamx2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmax(char scope, char topology, int m, int n, double* a,
  int lda, int rdest, int cdest) const
{ Cdgamx2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmax(char scope, int m, int n, double* a, int lda) const
{ dmax(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmax(int m, int n, double* a, int lda) const
{ dmax('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmin(char scope, char topology, int m, int n, double* a,
   int lda, int* ra, int* ca, int rcflag, int rdest, int cdest) const
{ Cdgamn2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmin(char scope, char topology, int m, int n, double* a,
  int lda, int rdest, int cdest) const
{ Cdgamn2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmin(char scope, int m, int n, double* a, int lda) const
{ dmin(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dmin(int m, int n, double* a, int lda) const
{ dmin('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_send(char scope, char topology, int m, int n,
  double* a,int lda) const
{ Cdgebs2d(ictxt_,&scope,&topology,m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_send(char scope, int m, int n, double* a, int lda) const
{ dbcast_send(scope,' ',m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_send(int m, int n, double* a, int lda) const
{ dbcast_send('A',' ',m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_recv(char scope, char topology, int m, int n,
  double* a, int lda, int rsrc, int csrc) const
{ Cdgebr2d(ictxt_,&scope,&topology,m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_recv(char scope,
  int m, int n, double* a, int lda, int rsrc, int csrc) const
{ dbcast_recv(scope,' ',m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::dbcast_recv(int m, int n, double* a, int lda,
  int rsrc,int csrc) const
{ dbcast_recv('A',' ',m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::isend(int m, int n, int* a, int lda, int rdest, int cdest)
  const { Cigesd2d(ictxt_,m,n,a,lda,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::irecv(int m, int n, int* a, int lda, int rsrc, int csrc) const
{ Cigerv2d(ictxt_,m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::isum(char scope, char topology, int m, int n, int* a, int lda,
  int rdest, int cdest) const
{ Cigsum2d(ictxt_,&scope,&topology,m,n,a,lda,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::isum(char scope, int m, int n, int* a, int lda) const
{ isum(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::isum(int m, int n, int* a, int lda) const
{ isum('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::imax(char scope, char topology, int m, int n, int* a, int lda,
  int* ra, int* ca, int rcflag, int rdest, int cdest) const
{ Cigamx2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::imax(char scope, char topology, int m, int n, int* a, int lda,
  int rdest, int cdest) const
{ Cigamx2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::imax(char scope, int m, int n, int* a, int lda) const
{ imax(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::imax(int m, int n, int* a, int lda) const
{ imax('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::imin(char scope, char topology, int m, int n, int* a, int lda,
  int* ra, int* ca, int rcflag, int rdest, int cdest) const
{ Cigamn2d(ictxt_,&scope,&topology,m,n,a,lda,ra,ca,rcflag,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::imin(char scope, char topology, int m, int n, int* a, int lda,
  int rdest, int cdest) const
{ Cigamn2d(ictxt_,&scope,&topology,m,n,a,lda,(int*)0,(int*)0,-1,rdest,cdest); }

////////////////////////////////////////////////////////////////////////////////
void Context::imin(char scope, int m, int n, int* a, int lda) const
{ imin(scope,' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::imin(int m, int n, int* a, int lda) const
{ imin('A',' ',m,n,a,lda,-1,-1); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_send(char scope, char topology,
  int m, int n, int* a,int lda) const
{ Cigebs2d(ictxt_,&scope,&topology,m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_send(char scope, int m, int n, int* a, int lda) const
{ ibcast_send(scope,' ',m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_send(int m, int n, int* a, int lda) const
{ ibcast_send('A',' ',m,n,a,lda); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_recv(char scope, char topology, int m, int n, int* a,
  int lda, int rsrc, int csrc) const
{ Cigebr2d(ictxt_,&scope,&topology,m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_recv(char scope,
  int m, int n, int* a, int lda, int rsrc, int csrc) const
{ ibcast_recv(scope,' ',m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::ibcast_recv(int m, int n, int* a, int lda,
  int rsrc,int csrc) const
{ ibcast_recv('A',' ',m,n,a,lda,rsrc,csrc); }

////////////////////////////////////////////////////////////////////////////////
void Context::string_send(std::string& s, int rdest, int cdest) const
{
  int len = s.size();
  isend(1,1,&len,1,rdest,cdest);
  int ilen = len/(sizeof(int)/sizeof(char));
  if ( len%(sizeof(int)/sizeof(char)) != 0 ) ilen++;
  int* ibuf = new int[ilen];
  s.copy((char*)ibuf,std::string::npos);
  isend(ilen,1,ibuf,ilen,rdest,cdest);
  delete [] ibuf;
}

////////////////////////////////////////////////////////////////////////////////
void Context::string_recv(std::string& s, int rsrc, int csrc) const
{
  int len = -1;
  irecv(1,1,&len,1,rsrc,csrc);
  int ilen = len/(sizeof(int)/sizeof(char));
  if ( len%(sizeof(int)/sizeof(char)) != 0 ) ilen++;
  int* ibuf = new int[ilen];
  irecv(ilen,1,ibuf,ilen,rsrc,csrc);
  s.resize(len);
  s.assign((char*)ibuf,len);
  delete [] ibuf;
}

////////////////////////////////////////////////////////////////////////////////
void Context::string_bcast(std::string& s, int isrc) const
{
  int len;
  if ( mype() == isrc )
  {
    len = s.length();
  }
  MPI_Bcast(&len,1,MPI_INT,isrc,comm());
  char* buf = new char[len+1];
  // s is initialized only on task isrc
  if ( mype() == isrc )
  {
    s.copy(buf,std::string::npos);
    buf[len]=0;
    assert(buf[len]=='\0');
  }
  MPI_Bcast(buf,len+1,MPI_CHAR,isrc,comm());
  s = buf;
  delete [] buf;
}

////////////////////////////////////////////////////////////////////////////////
void Context::print(ostream& os) const
{
  if ( active_ )
  {
    os << " " << nprow_ << "x" << npcol_ << " ";
    os << "{ ";
    for ( int ir = 0; ir < nprow_; ir++ )
    {
      os << "{ ";
      for ( int ic = 0; ic < npcol_; ic++ )
      {
        os << pmap_[ir+nprow_*ic] << " ";
      }
      os << "} ";
    }
    //os << "} (ictxt=" << ictxt_ << ")" << endl;
    os << "} (ictxt=" << ictxt_ << ")";
    int handle;
    Cblacs_get(ictxt_,10,&handle);
    os << "(handle=" << handle << ")";
    os << "(comm=" << comm_<< ")" << endl;
  }
  else
  {
    os << " (inactive)" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const Context& c)
{
  c.print(os);
  return os;
}
