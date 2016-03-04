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
// Context.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef CONTEXT_H
#define CONTEXT_H

#include <iosfwd>
#include <string>
#include <vector>
#include <cassert>
#include "blacs.h"
#include <mpi.h>

struct ContextRep
{
  private:

  int ictxt_;
  int myrow_;
  int mycol_;
  int nprow_;
  int npcol_;
  int size_;
  int myproc_;
  int mype_;
  bool onpe0_;
  bool active_;

  std::vector<int> pmap_;
  MPI_Comm comm_;

  // keep assignment and copy constructors private
  ContextRep& operator=(const ContextRep& c);
  ContextRep(const ContextRep& c);

  public:

  int ictxt() const { return ictxt_; }
  int myrow() const { return myrow_; }
  int mycol() const { return mycol_; }
  int nprow() const { return nprow_; }
  int npcol() const { return npcol_; }

  // number of processes in the context
  // returns -1 if current process is not part of this context
  int size() const { return size_; }
  // position of current process in row-major order
  // returns -1 if current process is not part of this context
  int myproc() const { return myproc_; }
  int mype() const { return mype_; }
  int pmap(int irow, int icol) const { return pmap_[irow+nprow_*icol]; }

  bool onpe0(void) const { return onpe0_; }
  bool active(void) const { return active_; }
  void abort(int ierr) const { Cblacs_abort(ictxt_,ierr); }
  void barrier(void) const { Cblacs_barrier(ictxt_,"A"); }
  void barrier(char scope) const { Cblacs_barrier(ictxt_,&scope); }

  void dsend(int m, int n, double* a, int lda, int rdest, int cdest) const;
  void drecv(int m, int n, double* a, int lda, int rsrc, int csrc) const;
  void dsum(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const;
  void dmax(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const;
  void dmax(char scope, char topology, int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmin(char scope, char topology, int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmin(char scope, char topology, int m, int n, double* a, int lda,
            int rdest, int cdest) const;
  void dbcast_send(char scope, char topology,
                   int m, int n, double* a,int lda) const;
  void dbcast_recv(char scope, char topology, int m, int n, double* a, int lda,
               int rsrc, int csrc) const;
  void isend(int m, int n, int* a, int lda, int rdest, int cdest) const;
  void irecv(int m, int n, int* a, int lda, int rsrc, int csrc) const;
  void isum(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const;
  void imax(char scope, char topology, int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imax(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const;
  void imin(char scope, char topology, int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imin(char scope, char topology, int m, int n, int* a, int lda,
            int rdest, int cdest) const;
  void ibcast_send(char scope, char topology,
                   int m, int n, int* a,int lda) const;
  void ibcast_recv(char scope, char topology, int m, int n, int* a, int lda,
               int rsrc, int csrc) const;
  void string_send(std::string& s, int rdest, int cdest) const;
  void string_recv(std::string& s, int rsrc, int csrc) const;
  void string_bcast(std::string& s, int isrc) const;

  bool operator==(const ContextRep& c) const
  { return (ictxt_==c.ictxt());}

  MPI_Comm comm(void) const { return comm_; }

  // Constructors

  // construct a single-row ContextRep
  explicit ContextRep(MPI_Comm comm);

  // global ContextRep of size nprow * npcol with column major order
  explicit ContextRep(MPI_Comm comm, int nprow, int npcol);

  ~ContextRep();

  void print(std::ostream& os) const;
};

class Context
{
  private:

  ContextRep* rep;
  int* pcount;

  public:

  int ictxt() const;
  int myrow() const;
  int mycol() const;
  int nprow() const;
  int npcol() const;

  // number of processes in the context
  // returns -1 if current process is not part of this context
  int size() const;
  // position of current process in row-major order
  // returns -1 if current process is not part of this context
  int myproc() const;
  int mype() const;
  // process id of process (irow,icol)
  int pmap(int irow, int icol) const;

  bool onpe0(void) const;
  bool active(void) const;
  operator bool() const { return active(); }
  void abort(int ierr) const;
  void barrier(void) const;
  void barrier(char scope) const;

  // double communications
  void dsend(int m, int n, double* a, int lda, int rdest, int cdest) const;
  void drecv(int m, int n, double* a, int lda, int rsrc, int csrc) const;

  void dsum(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dsum(char scope, int m, int n, double* a, int lda) const;
  void dsum(int m, int n, double* a, int lda) const;

  void dmax(char scope, char topology,
            int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmax(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dmax(char scope, int m, int n, double* a, int lda) const;
  void dmax(int m, int n, double* a, int lda) const;

  void dmin(char scope, char topology,
            int m, int n, double* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void dmin(char scope, char topology,
            int m, int n, double* a, int lda, int rdest, int cdest) const;
  void dmin(char scope, int m, int n, double* a, int lda) const;
  void dmin(int m, int n, double* a, int lda) const;

  void dbcast_send(char scope, char topology,
                   int m, int n, double* a, int lda) const;
  void dbcast_send(char scope, int m, int n, double* a, int lda) const;
  void dbcast_send(int m, int n, double* a, int lda) const;

  void dbcast_recv(char scope, char topology,
               int m, int n, double* a, int lda, int rsrc, int csrc) const;
  void dbcast_recv(char scope, int m, int n, double* a,
                   int lda,int rsrc, int csrc) const;
  void dbcast_recv(int m, int n, double* a, int lda,int rsrc, int csrc) const;

  // integer communications
  void isend(int m, int n, int* a, int lda, int rdest, int cdest) const;
  void irecv(int m, int n, int* a, int lda, int rsrc, int csrc) const;
  void isum(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void isum(char scope, int m, int n, int* a, int lda) const;
  void isum(int m, int n, int* a, int lda) const;

  void imax(char scope, char topology,
            int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imax(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void imax(char scope, int m, int n, int* a, int lda) const;
  void imax(int m, int n, int* a, int lda) const;

  void imin(char scope, char topology,
            int m, int n, int* a, int lda,
            int* ra, int* ca, int rcflag, int rdest, int cdest) const;
  void imin(char scope, char topology,
            int m, int n, int* a, int lda, int rdest, int cdest) const;
  void imin(char scope, int m, int n, int* a, int lda) const;
  void imin(int m, int n, int* a, int lda) const;

  void ibcast_send(char scope, char topology,
                   int m, int n, int* a, int lda) const;
  void ibcast_send(char scope, int m, int n, int* a, int lda) const;
  void ibcast_send(int m, int n, int* a, int lda) const;

  void ibcast_recv(char scope, char topology,
                   int m, int n, int* a, int lda, int rsrc, int csrc) const;
  void ibcast_recv(char scope, int m, int n,
                   int* a, int lda,int rsrc, int csrc) const;
  void ibcast_recv(int m, int n, int* a, int lda,int rsrc, int csrc) const;

  void string_send(std::string& s, int rdest, int cdest) const;
  void string_recv(std::string& s, int rsrc, int csrc) const;
  void string_bcast(std::string& s, int isrc) const;

  bool operator==(const Context& ctxt) const;

  // MPI communicator for this context. Returns MPI_COMM_NULL if
  // this process is not part of the context
  MPI_Comm comm(void) const;

  // Constructors

  // single-row Context
  explicit Context(MPI_Comm comm) : rep(new ContextRep(comm)),
    pcount(new int(1)) {}

  // nprow * npcol Context
  explicit Context(MPI_Comm comm, int nprow, int npcol):
    rep(new ContextRep(comm,nprow,npcol)), pcount(new int(1)) {}

//  Context(ContextRep* pp) : rep(pp), pcount(new int(1)) {}

  Context(const Context& c) : rep(c.rep), pcount(c.pcount) { (*pcount)++; }

  void print(std::ostream& os) const;

  Context& operator=(const Context& c)
  {
    if ( rep == c.rep ) return *this;
    if ( --(*pcount) == 0 )
    {
      delete rep;
      delete pcount;
    }
    rep = c.rep;
    pcount = c.pcount;
    (*pcount)++;
    return *this;
  }

  ~Context(void)
  {
    if ( pcount == 0 )
    {
      std::cerr << "~Context: pcount = 0\n";
    }
    if ( --(*pcount) == 0 )
    {
      delete rep;
      delete pcount;
    }
  }
};

std::ostream& operator << ( std::ostream& os, const Context& ctxt );

#endif
