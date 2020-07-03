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
// MPIdata.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MPIDATA_H
#define MPIDATA_H

#include <mpi.h>
class MPIdata
{
  static MPI_Comm comm_;
  static MPI_Comm g_comm_;
  static MPI_Comm st_comm_;
  static MPI_Comm sp_comm_;
  static MPI_Comm kp_comm_;
  static MPI_Comm sd_comm_;
  static int rank_;
  static int size_;
  static bool onpe0_;

  public:
  static const MPI_Comm& comm(void) { return comm_; }
  static const MPI_Comm& g_comm(void) { return g_comm_; }
  static const MPI_Comm& st_comm(void) { return st_comm_; }
  static const MPI_Comm& sp_comm(void) { return sp_comm_; }
  static const MPI_Comm& kp_comm(void) { return kp_comm_; }
  static const MPI_Comm& sd_comm(void) { return sd_comm_; }
  static int rank(void) { return rank_; }
  static int size(void) { return size_; }
  static bool onpe0(void) { return onpe0_; }
  static void set(int ngb, int nstb = 0, int nspb = 0, int nkpb = 0);
};
#endif
