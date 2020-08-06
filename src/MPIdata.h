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
  private:
  static MPI_Comm comm_;
  static MPI_Comm g_comm_;
  static MPI_Comm st_comm_;
  static MPI_Comm kp_comm_;
  static MPI_Comm sp_comm_;
  static MPI_Comm sd_comm_;
  static MPI_Comm kp_sp_comm_;
  static MPI_Comm st_kp_sp_comm_;
  static int rank_;
  static int size_;
  static bool onpe0_;

  static int ngb_;
  static int nstb_;
  static int nkpb_;
  static int nspb_;

  static int igb_;
  static int istb_;
  static int ikpb_;
  static int ispb_;

  static int sd_rank_;
  static int sd_size_;
  static int kp_sp_rank_;
  static int kp_sp_size_;
  static int st_kp_sp_rank_;
  static int st_kp_sp_size_;

  // private ctor
  MPIdata(void);

  public:
  static const MPI_Comm& comm(void) { return comm_; }
  static const MPI_Comm& g_comm(void) { return g_comm_; }
  static const MPI_Comm& st_comm(void) { return st_comm_; }
  static const MPI_Comm& kp_comm(void) { return kp_comm_; }
  static const MPI_Comm& sp_comm(void) { return sp_comm_; }
  static const MPI_Comm& sd_comm(void) { return sd_comm_; }
  static const MPI_Comm& kp_sp_comm(void) { return kp_sp_comm_; }
  static const MPI_Comm& st_kp_sp_comm(void) { return st_kp_sp_comm_; }

  static int rank(void) { return rank_; }
  static int size(void) { return size_; }
  static bool onpe0(void) { return onpe0_; }

  static int ngb(void)  { return ngb_; }
  static int nstb(void) { return nstb_; }
  static int nkpb(void) { return nkpb_; }
  static int nspb(void) { return nspb_; }

  static int igb(void)  { return igb_; }
  static int istb(void) { return istb_; }
  static int ikpb(void) { return ikpb_; }
  static int ispb(void) { return ispb_; }

  static int sd_rank(void) { return sd_rank_; }
  static int sd_size(void) { return sd_size_; }
  static int kp_sp_rank(void) { return kp_sp_rank_; }
  static int kp_sp_size(void) { return kp_sp_size_; }
  static int st_kp_sp_rank(void) { return st_kp_sp_rank_; }
  static int st_kp_sp_size(void) { return st_kp_sp_size_; }

  static void set(int ngb, int nstb = 1, int nkpb = 1, int nspb = 1);
};
#endif
