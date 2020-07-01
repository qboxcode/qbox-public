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
namespace MPIdata
{
  extern MPI_Comm comm;
  extern MPI_Comm g_comm;
  extern MPI_Comm st_comm;
  extern MPI_Comm sp_comm;
  extern MPI_Comm kp_comm;
  extern MPI_Comm sd_comm;
  extern int rank;
  extern int size;
  extern bool onpe0;
};
#endif
