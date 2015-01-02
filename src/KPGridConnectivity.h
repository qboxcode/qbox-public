////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2010 The Regents of the University of California
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
// KPGridConnectivity.h
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include "Sample.h"
#include "SlaterDet.h"

using namespace std;

#ifndef KPGRIDCONNECTIVITY_H
#define KPGRIDCONNECTIVITY_H

class KPConnectivity
{
  private:

  const Sample &s_;
  int DimX_;
  int DimY_;
  int DimZ_;
  int nkpoints_;
  int nStateMax_;
  int ConnectivityComplete_;
  double kdist_tol_;

  std::vector<int> first_neighbour_kx_;
  std::vector<int> first_neighbour_ky_;
  std::vector<int> first_neighbour_kz_;
  std::vector<int> second_neighbour_kx_;
  std::vector<int> second_neighbour_ky_;
  std::vector<int> second_neighbour_kz_;

  std::vector<int> first_symmetric_kx_;
  std::vector<int> first_symmetric_ky_;
  std::vector<int> first_symmetric_kz_;
  std::vector<int> second_symmetric_kx_;
  std::vector<int> second_symmetric_ky_;
  std::vector<int> second_symmetric_kz_;

  std::vector<int> local_ig_overlap_;
  std::vector<int> first_ig_overlap_kx_;
  std::vector<int> first_ig_overlap_ky_;
  std::vector<int> first_ig_overlap_kz_;
  std::vector<int> second_ig_overlap_kx_;
  std::vector<int> second_ig_overlap_ky_;
  std::vector<int> second_ig_overlap_kz_;

  std::vector<double> first_distance_kx_;
  std::vector<double> first_distance_ky_;
  std::vector<double> first_distance_kz_;
  std::vector<double> second_distance_kx_;
  std::vector<double> second_distance_ky_;
  std::vector<double> second_distance_kz_;

  std::vector<double> local_overlap_;
  std::vector<double> first_overlap_kx_;
  std::vector<double> first_overlap_ky_;
  std::vector<double> first_overlap_kz_;
  std::vector<double> second_overlap_kx_;
  std::vector<double> second_overlap_ky_;
  std::vector<double> second_overlap_kz_;

  std::vector<double> local_send_buff_;
  std::vector<double> first_send_buff_kx_;
  std::vector<double> first_send_buff_ky_;
  std::vector<double> first_send_buff_kz_;
  std::vector<double> second_send_buff_kx_;
  std::vector<double> second_send_buff_ky_;
  std::vector<double> second_send_buff_kz_;

  std::vector<double> integral_kx_;
  std::vector<double> integral_ky_;
  std::vector<double> integral_kz_;

  std::vector<double> weight_;

  std::vector<D3vector> first_T_overlap_kx_;
  std::vector<D3vector> first_T_overlap_ky_;
  std::vector<D3vector> first_T_overlap_kz_;
  std::vector<D3vector> second_T_overlap_kx_;
  std::vector<D3vector> second_T_overlap_ky_;
  std::vector<D3vector> second_T_overlap_kz_;

  MPI_Comm comm_;

  MPI_Request send_request_local_;
  MPI_Request send_request_first_kx_;
  MPI_Request send_request_first_ky_;
  MPI_Request send_request_first_kz_;
  MPI_Request send_request_second_kx_;
  MPI_Request send_request_second_ky_;
  MPI_Request send_request_second_kz_;

  MPI_Request recv_request_local_;
  MPI_Request recv_request_first_kx_;
  MPI_Request recv_request_first_ky_;
  MPI_Request recv_request_first_kz_;
  MPI_Request recv_request_second_kx_;
  MPI_Request recv_request_second_ky_;
  MPI_Request recv_request_second_kz_;

  public:

  // constructor
  KPConnectivity(const Sample& s_);

  // destructor
  ~KPConnectivity();

  // methods
  int Connection() { return ConnectivityComplete_; }
  void InitOverlaps();
  void SetOverlapIndices(Basis *vbasis_);
  void AddOverlap(int iKpi, int iKpj, int iLocStatei,
    complex<double> *valueDirect, complex<double> *valueSymmetric,
    double occupation);
  void cell_moved(void);
  void StartPermutation(int iKp, int iSendTo, int iRecvFr);
  void EndPermutation();

  double overlaps_local(int iKp, int iLocState);
  double overlaps_first_kx(int iKp, int iLocState);
  double overlaps_first_ky(int iKp, int iLocState);
  double overlaps_first_kz(int iKp, int iLocState);
  double overlaps_second_kx(int iKp, int iLocState);
  double overlaps_second_ky(int iKp, int iLocState);
  double overlaps_second_kz(int iKp, int iLocState);

  double distance_first_kx(int iKp);
  double distance_first_ky(int iKp);
  double distance_first_kz(int iKp);
  double distance_second_kx(int iKp);
  double distance_second_ky(int iKp);
  double distance_second_kz(int iKp);

  double integral_kx(int iKp);
  double integral_ky(int iKp);
  double integral_kz(int iKp);

  double weight(int iKp);
};
#endif
