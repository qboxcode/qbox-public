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

// global communicator
MPI_Comm MPIdata::comm;
int MPIdata::rank;     // global rank of this process
int MPIdata::size;     // global number of processes
bool MPIdata::onpe0;   // task 0 of global comm
// subcommunicators
MPI_Comm MPIdata::g_comm;   // G vector comm
MPI_Comm MPIdata::st_comm;  // states comm
MPI_Comm MPIdata::sp_comm;  // spin comm
MPI_Comm MPIdata::kp_comm;  // kpoint comm
MPI_Comm MPIdata::sd_comm;  // Slater det comm
