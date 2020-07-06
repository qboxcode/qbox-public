////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2020 The Regents of the University of California
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
// MPIdata.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include<iostream>
#include<cassert>
using namespace std;

// global communicator
MPI_Comm MPIdata::comm_;
// subcommunicators
MPI_Comm MPIdata::g_comm_;   // G vector comm
MPI_Comm MPIdata::st_comm_;  // states comm
MPI_Comm MPIdata::sp_comm_;  // spin comm
MPI_Comm MPIdata::kp_comm_;  // kpoint comm
MPI_Comm MPIdata::sd_comm_;  // Slater det comm
int MPIdata::rank_;     // global rank of this process
int MPIdata::size_;     // global number of processes
bool MPIdata::onpe0_;   // task 0 of global comm

void MPIdata::set(int ngb, int nstb, int nspb, int nkpb)
{
  MPI_Comm_size(MPI_COMM_WORLD,&size_);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
  onpe0_ = ( rank_ == 0 );

  assert(ngb*nstb*nspb*nkpb == size_);

  // check that all numbers of blocks are positive
  assert(ngb>0);
  assert(nstb>0);
  assert(nspb>0);
  assert(nkpb>0);

  int ndims=4;
  // Note: order of dimensions for Cart_comm to have contiguous tasks
  // along the ngb dimension
  int dims[4] = { nkpb, nspb, nstb, ngb };
  // plan for cyclic rotation of states blocks
  int periods[4] = { 0, 0, 1, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm_);

  // create subcommunicators
  // G vector communicator
  int g_remain_dims[4] = { 0, 0, 0, 1 };
  MPI_Cart_sub(comm_,g_remain_dims,&g_comm_);

  // states communicator
  int st_remain_dims[4] = { 0, 0, 1, 0 };
  MPI_Cart_sub(comm_,st_remain_dims,&st_comm_);

  // spin communicator
  int sp_remain_dims[4] = { 0, 1, 0, 0 };
  MPI_Cart_sub(comm_,sp_remain_dims,&sp_comm_);

  // kpoint communicator
  int kp_remain_dims[4] = { 1, 0, 0, 0 };
  MPI_Cart_sub(comm_,kp_remain_dims,&kp_comm_);

  // Slater determinant communicator
  int sd_remain_dims[4] = { 0, 0, 1, 1 };
  MPI_Cart_sub(comm_,sd_remain_dims,&sd_comm_);
}
