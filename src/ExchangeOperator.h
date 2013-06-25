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
// ExchangeOperator.h
//
////////////////////////////////////////////////////////////////////////////////

#include "Sample.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "KPGridConnectivity.h"

#ifndef EXCHANGEOPERATOR_H
#define EXCHANGEOPERATOR_H

class Bisection;
class ExchangeOperator
{
  private:

  Sample& s_;
  // exchange energy
  double eex_;
  // constant of support function for exchange integration
  double rcut_;
  // mixing coefficient for exchange energy and dwf accumulation
  double HFCoeff_;

  // HF stress tensor
  std::valarray<double> sigma_exhf_;

  // reference wf and dwf for non scf iteration
  Wavefunction wf0_;
  Wavefunction dwf0_;
  // copy of wf
  Wavefunction wfc_;

  double compute_exchange_for_general_case_(Sample* s, Wavefunction* dwf,
    bool compute_stress);
  double compute_exchange_at_gamma_(const Wavefunction &wf, Wavefunction* dwf,
    bool compute_stress);
  void   apply_VXC_(double mix, Wavefunction& wf_ref,
    Wavefunction& dwf_ref, Wavefunction& dwf);

  // Connectivity of the kpoint Grid
  KPConnectivity KPGridPerm_;
  KPConnectivity KPGridStat_;

  // basis for pair densities
  Basis* vbasis_;
  int np0v_, np1v_, np2v_;
  bool gamma_only_;

  // Fourier transform for pair densities
  FourierTransform* vft_;
  int np012loc_;

  // Fourier transform for states and forces at gamma
  FourierTransform*  wft_;
  FourierTransform* dwft_;

  // pair densities
  vector<complex<double> > rhog1_;
  vector<complex<double> > rhog2_;
  valarray<complex<double> > rhor1_;
  valarray<complex<double> > rhor2_;

  // G vectors
  valarray<double> qpG21_;
  valarray<double> qpG22_;
  valarray<double> qpG2i1_;
  valarray<double> qpG2i2_;
  valarray<double> G2_;
  valarray<double> G2i_;

  // numbers of states
  int nLocalStates_;
  int nMaxLocalStates_;

  // exchange energies
  vector<double> exchange_ki_;
  vector<double> exchange_kj_;
  vector<double> send_buf_exchange_;
  vector<double> send_buf_occupation_;
  vector<vector<double> > exchange_;

  // states arrays
  valarray<complex<double> > tmp_;
  vector<complex<double> > state_kpi_;
  vector<complex<double> > send_buf_states_;
  vector<valarray<complex<double> > > statei_;
  vector<valarray<complex<double> > > statej_;

  // occupation numbers
  vector<double> occ_ki_;
  vector<double> occ_kj_;

  // forces array
  vector<complex<double> > force_kpi_;
  vector<complex<double> > send_buf_forces_;
  vector<valarray<complex<double> > > dstatei_;
  vector<valarray<complex<double> > > dstatej_;

  // buffer
  vector<complex<double> > buffer_dstate_;
  vector<complex<double> > buffer_forces_1_;
  vector<complex<double> > buffer_forces_2_;

  // contexts and communicators
  Context gcontext_;
  MPI_Comm vcomm_;

  // Communications
  int colSendTo_;
  int colRecvFr_;
  int iSendTo_;
  int iRecvFr_;

  int nStatesKpi_;
  int nNextStatesKpi_;

  // communicator
  MPI_Comm comm_;

  // MPI communications request
  MPI_Request send_request_NumberOfStates_;
  MPI_Request recv_request_NumberOfStates_;
  MPI_Request send_request_States_;
  MPI_Request recv_request_States_;
  MPI_Request send_request_Forces_;
  MPI_Request recv_request_Forces_;
  MPI_Request send_request_Exchange_;
  MPI_Request recv_request_Exchange_;
  MPI_Request send_request_Occupation_;
  MPI_Request recv_request_Occupation_;

  // send/reception flags are used in case of 0 states permutation
  int wait_send_states_;
  int wait_recv_states_;
  int wait_send_forces_;
  int wait_recv_forces_;
  int wait_send_energies_;
  int wait_recv_energies_;
  int wait_send_occupations_;
  int wait_recv_occupations_;

  // communication functions
  void InitPermutation(void);
  void FreePermutation(void);
  void SetNextPermutationStateNumber(void);
  void StartStatesPermutation(int mloc);
  void CompleteReceivingStates(int iRotationStep);
  void CompleteSendingStates(int iRotationStep);
  void StartForcesPermutation(int mloc);
  void CompleteReceivingForces(int iRotationStep);
  void CompleteSendingForces(int iRotationStep);
  void StartEnergiesPermutation(void);
  void CompleteReceivingEnergies(int iRotationStep);
  void CompleteSendingEnergies(int iRotationStep);
  void StartOccupationsPermutation(void);
  void CompleteReceivingOccupations(int iRotationStep);
  void CompleteSendingOccupations(int iRotationStep);

  // bisection
  bool use_bisection_;
  vector<Bisection*> bisection_;
  vector<DoubleMatrix*> uc_;
  vector<long int> localization_;

  public:

  // constructor
  ExchangeOperator(Sample& s_, double HFCoeff);

  // destructor
  ~ExchangeOperator();

  // parameters
  void setmixCoeff(double value) { HFCoeff_ = value; };
  double HFCoeff() { return HFCoeff_; };

  // exchange energy and forces computation
  double eex() { return eex_; };
  double update_energy(bool compute_stress);
  double update_operator(bool compute_stress);
  double apply_operator(Wavefunction& dwf);
  void add_stress (std::valarray<double> & sigma_exc);
};

class ExchangeOperatorException
{
  public:
  std::string msg;
  ExchangeOperatorException(std::string s) : msg(s) {}
};
#endif
