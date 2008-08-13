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
// SharedFilePtr.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SharedFilePtr.h,v 1.2 2008-08-13 06:39:43 fgygi Exp $

#ifndef SHAREDFILEPTR_H
#define SHAREDFILEPTR_H

#include "mpi.h"

class SharedFilePtr
{
  private:

  MPI_Comm comm_;
  MPI_File& fh_;
  long long int offset_;

  public:

  MPI_File& file(void) { return fh_; }
  long long int offset(void) const { return offset_; }
  MPI_Offset mpi_offset(void) const { return (MPI_Offset) offset_; }
  void sync(void)
  {
    // set all offsets to the largest offset
    long long int s_off = offset_;
    MPI_Allreduce(&s_off,&offset_,1,MPI_LONG_LONG,MPI_MAX,comm_);
  }
  void set_offset(long long int off)
  {
    offset_ = off;
  }
  void advance(long long int dist)
  {
    offset_ += dist;
  }

  SharedFilePtr(MPI_Comm comm, MPI_File& fh) : comm_(comm),
                fh_(fh), offset_(0) {}
  ~SharedFilePtr() {}
};
#endif
