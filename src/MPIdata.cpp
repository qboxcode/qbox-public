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
MPI_Comm MPIdata::kp_comm_;  // kpoint comm
MPI_Comm MPIdata::sp_comm_;  // spin comm
MPI_Comm MPIdata::sd_comm_;  // Slater det comm
MPI_Comm MPIdata::kp_sp_comm_;  // kpoints+spin comm
MPI_Comm MPIdata::st_kp_sp_comm_;  // states+kpoints+spin comm
int MPIdata::rank_;     // global rank of this process
int MPIdata::size_;     // global number of processes
bool MPIdata::onpe0_;   // task 0 of global comm

int MPIdata::ngb_;      // size of g_comm
int MPIdata::nstb_;     // size of st_comm
int MPIdata::nkpb_;     // size of kp_comm
int MPIdata::nspb_;     // size of sp_comm

int MPIdata::igb_;      // rank in g_comm
int MPIdata::istb_;     // rank in st_comm
int MPIdata::ikpb_;     // rank in kp_comm
int MPIdata::ispb_;     // rank in sp_comm

int MPIdata::sd_rank_;  // rank in sd_comm
int MPIdata::sd_size_;  // size of sd_comm
int MPIdata::kp_sp_rank_;  // rank in kp_sp_comm
int MPIdata::kp_sp_size_;  // size of kp_sp_comm
int MPIdata::st_kp_sp_rank_;  // rank in st_kp_sp_comm
int MPIdata::st_kp_sp_size_;  // size of st_kp_sp_comm

void MPIdata::set(int ngb, int nstb, int nkpb, int nspb)
{
  MPI_Comm_size(MPI_COMM_WORLD,&size_);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
  onpe0_ = ( rank_ == 0 );

  assert(ngb*nstb*nkpb*nspb == size_);

  // check that all numbers of blocks are positive
  assert(ngb>0);
  assert(nstb>0);
  assert(nkpb>0);
  assert(nspb>0);

  int ndims=4;
  // Note: order of dimensions for Cart_comm to have contiguous tasks
  // along the ngb dimension
  int dims[4] = { nspb, nkpb, nstb, ngb };
  // plan for cyclic rotation of states blocks
  int periods[4] = { 0, 0, 1, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm_);

  // create subcommunicators
  // G vector communicator
  int g_remain_dims[4] = { 0, 0, 0, 1 };
  MPI_Cart_sub(comm_,g_remain_dims,&g_comm_);
  MPI_Comm_size(g_comm_,&ngb_);
  MPI_Comm_rank(g_comm_,&igb_);

  // states communicator
  int st_remain_dims[4] = { 0, 0, 1, 0 };
  MPI_Cart_sub(comm_,st_remain_dims,&st_comm_);
  MPI_Comm_size(st_comm_,&nstb_);
  MPI_Comm_rank(st_comm_,&istb_);

  // kpoint communicator
  int kp_remain_dims[4] = { 0, 1, 0, 0 };
  MPI_Cart_sub(comm_,kp_remain_dims,&kp_comm_);
  MPI_Comm_size(kp_comm_,&nkpb_);
  MPI_Comm_rank(kp_comm_,&ikpb_);

  // spin communicator
  int sp_remain_dims[4] = { 1, 0, 0, 0 };
  MPI_Cart_sub(comm_,sp_remain_dims,&sp_comm_);
  MPI_Comm_size(sp_comm_,&nspb_);
  MPI_Comm_rank(sp_comm_,&ispb_);

  // Slater determinant communicator
  int sd_remain_dims[4] = { 0, 0, 1, 1 };
  MPI_Cart_sub(comm_,sd_remain_dims,&sd_comm_);
  MPI_Comm_size(sd_comm_,&sd_size_);
  MPI_Comm_rank(sd_comm_,&sd_rank_);

  // kp_sp communicator
  int kp_sp_remain_dims[4] = { 1, 1, 0, 0 };
  MPI_Cart_sub(comm_,kp_sp_remain_dims,&kp_sp_comm_);
  MPI_Comm_size(kp_sp_comm_,&kp_sp_size_);
  MPI_Comm_rank(kp_sp_comm_,&kp_sp_rank_);

  // st_kp_sp communicator
  int st_kp_sp_remain_dims[4] = { 1, 1, 1, 0 };
  MPI_Cart_sub(comm_,st_kp_sp_remain_dims,&st_kp_sp_comm_);
  MPI_Comm_size(st_kp_sp_comm_,&st_kp_sp_size_);
  MPI_Comm_rank(st_kp_sp_comm_,&st_kp_sp_rank_);
}
