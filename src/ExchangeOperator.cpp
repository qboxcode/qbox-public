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
// ExchangeOperator.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// Screened exchange operator
// Diagonal screening is used with interaction potential
// vint(r) = alpha_sx  * erf(mu_sx*r)/r + beta_sx * erfc(mu_sx*r)/r
// The parameters alpha_sx, beta_sx and mu_sx are provided to the constructor
//
// Special cases:
// alpha_sx = beta_sx = 0: no interaction
// alpha_sx = beta_sx: Coulomb potential with prefactor alpha_sx (= beta_sx)
// alpha_sx = 0, beta_sx = 0.25, mu = 0.11: HSE
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <iomanip>
#include <bitset>
#include <algorithm>
#include "VectorLess.h"

#include "ExchangeOperator.h"
#include "Bisection.h"
#include "MPIdata.h"
#include "cout0.h"

using namespace std;

enum MessageType { NumberOfStates, Occupation, Exchange, Forces, States };

////////////////////////////////////////////////////////////////////////////////
ExchangeOperator::ExchangeOperator( Sample& s, double alpha_sx,
  double beta_sx, double mu_sx ) :
  s_(s), wf0_(s.wf), dwf0_(s.wf), wfc_(s.wf),
  alpha_sx_(alpha_sx), beta_sx_(beta_sx), mu_sx_(mu_sx), wft_(0)
{
  // check validity of the values of alpha_sx, beta_sx, mu_sx
  // if alpha_sx == beta_sx (Coulomb potential) mu_sx is not used
  // if alpha_sx != beta_sx, mu_sx_ must be positive
  if ( alpha_sx != beta_sx ) assert( mu_sx > 0.0 );

  eex_ = 0.0; // exchange energy
  rcut_ = 1.0;  // constant of support function for exchange integration

  sigma_exhf_.resize(6);

  // check if the only kpoint is the gamma point:
  gamma_only_ = ( s_.wf.nkp()==1 ) && ( norm2(s_.wf.kpoint(0)) == 0.0 );

  if ( gamma_only_ )
  {
    // create a real basis for the pair densities
    vbasis_ = new Basis(MPIdata::g_comm(), D3vector(0.0,0.0,0.0));
  }
  else
  {
    // create a complex basis
    //!! should avoid the finite k trick to get a complex basis at gamma
    vbasis_ = new Basis(MPIdata::g_comm(),
      D3vector(0.00000001,0.00000001,0.00000001));
  }
  // the size of the basis for the pair density should be
  // twice the size of the wave function basis
  vbasis_->resize(s_.wf.cell(),s_.wf.refcell(),4.0*s_.wf.ecut());

  // set the size for the r grid to be a product of small primes
  np0v_ = vbasis_->np(0)+2;
  np1v_ = vbasis_->np(1)+2;
  np2v_ = vbasis_->np(2)+2;
  while (!vbasis_->factorizable(np0v_)) np0v_ += 2;
  while (!vbasis_->factorizable(np1v_)) np1v_ += 2;
  while (!vbasis_->factorizable(np2v_)) np2v_ += 2;

  const int ngloc = vbasis_->localsize();
  // create Fourier transform object for densities
  vft_ = new FourierTransform(*vbasis_,np0v_,np1v_,np2v_);

  np012loc_ = vft_->np012loc();

  // allocate memory for densities in G space
  rhog1_.resize(ngloc);
  rhog2_.resize(ngloc);

  // allocate memory for densities in r space
  rhor1_.resize(np012loc_);
  rhor2_.resize(np012loc_);

  // if not only at gamma, allocate arrays q+G
  if ( !gamma_only_ )
  {
    int_pot1_.resize(ngloc);
    int_pot2_.resize(ngloc);
  }

  // get maximum number of states on a task
  if ( s_.wf.nspin()==1 )
  {
    const int isp_loc = s_.wf.isp_local(0);
    int nLocalStates_ = 0;
    if ( isp_loc >= 0 )
      nLocalStates_ = s_.wf.sd(isp_loc,0)->nstloc();
    MPI_Allreduce(&nLocalStates_,&nMaxLocalStates_,1,
      MPI_INT,MPI_MAX,MPIdata::comm());
  }
  else
  {
    int nLocalStates_ = 0;
    for ( int ispin = 0; ispin < 2; ++ispin )
    {
      const int isp_loc = s_.wf.isp_local(ispin);
      if ( isp_loc >= 0 )
        nLocalStates_ = max(nLocalStates_,s_.wf.sd(isp_loc,0)->nstloc());
    }
    MPI_Allreduce(&nLocalStates_,&nMaxLocalStates_,1,
      MPI_INT,MPI_MAX,MPIdata::comm());
  }

  // allocate memory for exchange energies
  exchange_ki_.resize(nMaxLocalStates_);
  exchange_kj_.resize(nMaxLocalStates_);
  send_buf_exchange_.resize(nMaxLocalStates_);
  send_buf_occupation_.resize(nMaxLocalStates_);

  // allocate memory for exchange
  exchange_.resize(s_.wf.nkp());
  for ( int iKp=0; iKp<s_.wf.nkp(); iKp++ )
  {
    // allocate memory for exchange energies of states of this kpoint.
    exchange_[iKp].resize(nMaxLocalStates_);
  }
  // get maximum number of g coeff per states
  int mlocMax = 0;
  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
    {
      SlaterDet* sd = s_.wf.sd(isp_loc,ikp_loc);
      ComplexMatrix& c = sd->c();
      if ( mlocMax < c.mloc() ) mlocMax = c.mloc();
    }
  }
  int mlocMax_tmp = mlocMax;
  MPI_Allreduce(&mlocMax_tmp,&mlocMax,1,MPI_INT,MPI_MAX,MPIdata::comm());

  // allocate memory for the copy of states of kpoint iKpi
  state_kpi_.resize( nMaxLocalStates_ * mlocMax );
  send_buf_states_.resize( nMaxLocalStates_ * mlocMax );

  // allocate buffers ( different if only at gamma or not )
  if ( gamma_only_ )
  {
    buffer_forces_1_.resize( mlocMax );
    buffer_forces_2_.resize( mlocMax );
  }
  else
  {
    buffer_dstate_.resize( mlocMax );
  }

  // allocate memory for the r coordinate expression
  // of each state of kpi and kpj
  statei_.resize(nMaxLocalStates_);
  statej_.resize(nMaxLocalStates_);
  for ( int i = 0; i < nMaxLocalStates_; i++ )
  {
    statei_[i].resize(np012loc_);
    statej_[i].resize(np012loc_);
  }

  use_bisection_ = s.ctrl.btHF > 0.0;

  // if only at gamma
  if ( gamma_only_ )
  {
    tmp_.resize(np012loc_);
    // allocate bisection object
    if ( use_bisection_ )
    {
      bisection_.resize(s_.wf.nsp_loc());
      uc_.resize(s_.wf.nsp_loc());
      for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
      {
        bisection_[isp_loc] = new Bisection(*s_.wf.sd(isp_loc,0),s_.ctrl.blHF);
        const ComplexMatrix& c = s_.wf.sd(isp_loc,0)->c();
        uc_[isp_loc] = new DoubleMatrix(c.context(),c.n(),c.n(),c.nb(),c.nb());
      }
    }
  }

  // allocate memory for occupation numbers of kpoint iKpi
  occ_ki_.resize(nMaxLocalStates_);
  occ_kj_.resize(nMaxLocalStates_);

  // allocate memory for the real space expression of the forces
  dstatei_.resize(nMaxLocalStates_);
  dstatej_.resize(nMaxLocalStates_);
  for ( int i = 0; i < nMaxLocalStates_; i++ )
  {
    dstatei_[i].resize(np012loc_);
    dstatej_[i].resize(np012loc_);
  }

  // allocate memory for the copy of forces of kpoint iKpi
  {
    force_kpi_.resize( nMaxLocalStates_ * mlocMax );
    send_buf_forces_.resize( nMaxLocalStates_ * mlocMax );
  }
}

////////////////////////////////////////////////////////////////////////////////
ExchangeOperator::~ExchangeOperator()
{
  // delete Fourier transform objects on states and forces
  delete wft_;
  // delete Fourier transform and basis for pair densities
  delete vft_;
  delete vbasis_;
  if ( use_bisection_ )
    for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
    {
      delete bisection_[isp_loc];
      delete uc_[isp_loc];
    }
}

////////////////////////////////////////////////////////////////////////////////
double ExchangeOperator::update_energy(bool compute_stress)
{
  if ( gamma_only_ )
    return eex_ = compute_exchange_at_gamma_(s_.wf, 0, compute_stress);
  else
    return eex_ = compute_exchange_for_general_case_(s_.wf, 0, compute_stress);
}

////////////////////////////////////////////////////////////////////////////////
double ExchangeOperator::update_operator(bool compute_stress)
{
  dwf0_.clear();

  // compute exchange energy and derivatives
  if ( gamma_only_ )
    eex_ = compute_exchange_at_gamma_(s_.wf, &dwf0_, compute_stress);
  else
    eex_ = compute_exchange_for_general_case_(s_.wf, &dwf0_, compute_stress);

  // wf0_ is kept as a reference state
  wf0_ = s_.wf;

  // return exchange energy
  return eex_;
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::apply_VXC_(double mix, Wavefunction& wf_ref,
  Wavefunction& dwf_ref, Wavefunction& dwf)
{
  // dwf += mix * ( |dwf_ref> <wf_ref|s_.wf> + |wf_ref><dwf_ref|s_.wf>
  // - |wf_ref><wf_ref|sigma_HF|wf_ref><wf_ref|s_.wf> )
  for ( int isp_loc = 0; isp_loc < s_.wf.nsp_loc(); ++isp_loc )
  {
    const int ispin = s_.wf.isp_global(isp_loc);
    const int nst = s_.wf.nst(ispin);
    for ( int ikp_loc = 0; ikp_loc < s_.wf.nkp_loc(); ++ikp_loc )
    {
      const Context &ctxt = s_.wf.sd(isp_loc,ikp_loc)->c().context();
      if ( s_.wf.sd(isp_loc,ikp_loc)->basis().real() )
      {
        DoubleMatrix c_proxy(s_.wf.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix dc_proxy(dwf.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix cref_proxy(wf_ref.sd(isp_loc,ikp_loc)->c());
        DoubleMatrix dcref_proxy(dwf_ref.sd(isp_loc,ikp_loc)->c());

        int nb = c_proxy.nb();
        DoubleMatrix matproj1(ctxt,nst,nst,nb,nb);
        DoubleMatrix matproj2(ctxt,nst,nst,nb,nb);
        DoubleMatrix matenergy(ctxt,nst,nst,nb,nb);

        // matproj1 = <wf_ref|wf> => matproj1
        matproj1.gemm('t','n',2.0,cref_proxy,c_proxy,0.0);
        matproj1.ger(-1.0,cref_proxy,0,c_proxy,0);

        // dwf += mix * |dwf_ref> * matproj1
        dc_proxy.gemm('n','n',mix,dcref_proxy,matproj1,1.0);

        // matenergy = <dpsi_ref|psi_ref>
        matenergy.gemm('t','n',2.0,dcref_proxy,cref_proxy,0.0);
        matenergy.ger(-1.0,dcref_proxy,0,cref_proxy,0);

        // matproj2 = - matenergy * matproj1
        matproj2.gemm('n','n',-1.0,matenergy,matproj1,0.0);

        // matproj2 += <dwf_ref|wf>
        matproj2.gemm('t','n',2.0,dcref_proxy,c_proxy,1.0);
        matproj2.ger(-1.0,dcref_proxy,0,c_proxy,0);

        // |dwf> += mix * |wf_ref> * matproj2
        dc_proxy.gemm('n','n',mix,cref_proxy,matproj2,1.0);
      }
      else // complex wave functions
      {
        ComplexMatrix &c(s_.wf.sd(isp_loc,ikp_loc)->c());
        ComplexMatrix &dc(dwf.sd(isp_loc,ikp_loc)->c());
        ComplexMatrix &cref(wf_ref.sd(isp_loc,ikp_loc)->c());
        ComplexMatrix &dcref(dwf_ref.sd(isp_loc,ikp_loc)->c());

        int nb = c.nb();
        ComplexMatrix matproj1(ctxt,nst,nst,nb,nb);
        ComplexMatrix matproj2(ctxt,nst,nst,nb,nb);
        ComplexMatrix matenergy(ctxt,nst,nst,nb,nb);

        // matproj1 = <wf_ref|wf>
        matproj1.gemm('c','n',1.0,cref,c,0.0);

        // dwf += mix * |dwf_ref> * matproj1
        dc.gemm('n','n',mix,dcref,matproj1,1.0);

        // matenergy = <dwf_ref|wf_ref>
        matenergy.gemm('c','n',1.0,dcref,cref,0.0);

        // matproj2 = - matenergy * matproj1
        matproj2.gemm('n','n',-1.0,matenergy,matproj1,0.0);

        // matproj2 += <dwf_ref|wf>
        matproj2.gemm('c','n',1.0,dcref,c,1.0);

        // |dpsi> += mix * |psi_ref> * matproj2
        dc.gemm('n','n',mix,cref,matproj2,1.0);
      }
    } // ikp_loc
  } // isp_loc
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::apply_operator(Wavefunction& dwf)
{
  // apply sigmaHF to s_.wf and store result in dwf
  // use the reference function wf0_ and reference sigma(wf) dwf0_
  apply_VXC_(1.0, wf0_, dwf0_, dwf);
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::add_stress(valarray<double>& sigma_exc)
{
  // add current value of the HF stress tensor to sigma_exc
  sigma_exc += sigma_exhf_;
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::cell_moved(void)
{
  vbasis_->resize( s_.wf.cell(),s_.wf.refcell(),4.0*s_.wf.ecut());
}

////////////////////////////////////////////////////////////////////////////////
// Exchange functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
double ExchangeOperator::compute_exchange_for_general_case_
  (const Wavefunction& wf, Wavefunction* dwf, bool compute_stress)
{

  cout << "ExchangeOperator:: general kpoints not implemented" << endl;
  assert(false);

  Timer tm;
  tm.start();

  const double omega = wf.cell().volume();
  const int nkpoints = wf.nkp();
  const int nspin = wf.nspin();

  // determine the spin factor for the pair densities:
  // 0.5 if 1 spin, 1 if nspin==2
  const double spinFactor = 0.5 * nspin;
  const double exfac = - ( 4.0 * M_PI / omega ) * spinFactor;

  // initialize total exchange energy
  double exchange_sum = 0.0;

#if 0
  // initialize stress
  sigma_exhf_ = 0.0;

  // loop on spins
  for ( int ispin=0; ispin < nspin; ispin++ )
  {
    vector<double> num_corr(nkpoints);
    vector<vector<double> > sigma_num_corr(nkpoints);
    for ( int iKpi = 0; iKpi < nkpoints; iKpi++ )
    {
      num_corr[iKpi]=0.0;
      sigma_num_corr[iKpi].resize(6,0.0);
    }

    for ( int iKpi = 0; iKpi < nkpoints; iKpi++ )
    {
      SlaterDet& sdi = *(wf.sd(ispin,iKpi));
      ComplexMatrix& ci = sdi.c();
      FourierTransform* wfti_ = new FourierTransform(sdi.basis(),
        np0v_,np1v_,np2v_);
      SlaterDet& dsdi = *(dwf->sd(ispin,iKpi));
      ComplexMatrix& dci = dsdi.c();

      nStatesKpi_=sdi.nstloc();

      // occupation numbers for kpoint i
      const double* occ = sdi.occ_ptr();
      for ( int i = 0; i<nStatesKpi_; i++ )
        occ_ki_[i] = occ[ci.jglobal(i)];

      // copy of the local states at kpoint iKpi
      const complex<double> *p = ci.cvalptr(0);
      for ( int i = 0; i < nStatesKpi_ * ci.mloc(); i++ )
        state_kpi_[i]=p[i];

      if (dwf)
      {
        for ( int i = 0; i < nStatesKpi_ * ci.mloc(); i++ )
          force_kpi_[i] = 0.0;
      }

      // initialize communications for the permutations
      InitPermutation();

      //  Start rotation of the states of kpoint i from this point
      for ( int iRotationStep=0; iRotationStep<sdi.context().npcol();
            iRotationStep++ )
      {
        CompleteReceivingStates(iRotationStep);

        // compute the r coordinate expression of each state of kpi
        for ( int i = 0; i < nStatesKpi_; i++ )
        {
          wfti_->backward( &state_kpi_[i*ci.mloc()], &(statei_[i])[0] );
        }

        CompleteSendingStates(iRotationStep);

        for ( int i = 0; i < nStatesKpi_ * ci.mloc(); i++ )
          send_buf_states_[i]=state_kpi_[i];

        if (dwf)
        {
          for ( int i = 0; i < nStatesKpi_; i++ )
            for ( int j = 0; j < np012loc_; j++ )
              dstatei_[i][j]=0.0;
        }

        // set number of states of next permutation step
        // (contained in nNextStatesKpi_)
        SetNextPermutationStateNumber();

        // start states permutation
        StartStatesPermutation(ci.mloc());

        CompleteReceivingOccupations(iRotationStep);

        // second loop over kpoints
        for ( int iKpj = iKpi; iKpj < nkpoints; iKpj++ )
        {
          SlaterDet& sdj = *(wf.sd(ispin,iKpj));
          FourierTransform* wftj_ = new FourierTransform(sdj.basis(),
            np0v_,np1v_,np2v_);
          ComplexMatrix& cj = sdj.c();

          for ( int i = 0; i < sdj.nstloc(); i++ )
            wftj_->backward(cj.cvalptr(i*cj.mloc()),&(statej_[i])[0]);

          SlaterDet& dsdj = *(dwf->sd(ispin,iKpj));
          ComplexMatrix& dcj = dsdj.c();

          if (dwf)
          {
            for ( int i = 0; i < dsdj.nstloc(); i++ )
              for ( int j = 0; j < np012loc_; j++ )
                dstatej_[i][j]=0.0;
          }

          // compute the differences dk between kpoints i and j
          // and their symmetric equivalent
          // dk1 = kpi-kpj
          D3vector dk1 =   wf.kpoint(iKpi) - wf.kpoint(iKpj);

          // set dk1 in absolute reciprocal coordinates to get q1
          D3vector q1 = dk1.x*sdi.basis().cell().b(0)
                     +  dk1.y*sdi.basis().cell().b(1)
                     +  dk1.z*sdi.basis().cell().b(2);

          // dk2 = kpi+kpj
          D3vector dk2 =   wf.kpoint(iKpi) + wf.kpoint(iKpj);

          // set dk2 in absolute reciprocal coordinates to get q2
          D3vector q2 = dk2.x*sdi.basis().cell().b(0)
                     +  dk2.y*sdi.basis().cell().b(1)
                     +  dk2.z*sdi.basis().cell().b(2);

          // compute correction term exp(-rcut_^2*(q+G)^2)/(q+G)^2
          // compute and store vint(q1+G) and vint(q2+G)
          double SumExpQpG2 = 0.0;
          double sigma_sumexp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

          const int ngloc = vbasis_->localsize();
          for ( int ig = 0; ig < ngloc; ig++ )
          {
            D3vector G(vbasis_->gx(ig+ngloc*0),
                       vbasis_->gx(ig+ngloc*1),
                       vbasis_->gx(ig+ngloc*2));
            D3vector q1pG(q1+G);
            D3vector q2pG(q2+G);
            const double q1pG2 = norm2(q1pG);
            int_pot1_[ig] = vint(q1pG2);

            const double q2pG2 = norm2(q2pG);
            int_pot2_[ig] = vint(q2pG2);

            // if iKpi=0 (first k point)
            // correction term: sum_(G,q) exp(-rcut_^2*|q+G|^2)/(q+G)^2
            if ( (iRotationStep==0) && (alpha_sx_ != 0.0))
            {
              const double rc2 = rcut_*rcut_;
              const double q1pG2i = ( q1pG2 > 0.0 ) ? 1.0 / q1pG2 : 0.0;
              const double q2pG2i = ( q2pG2 > 0.0 ) ? 1.0 / q2pG2 : 0.0;
              const double t1 = alpha_sx_ * (exp(-rc2*q1pG2) * q1pG2i );
              const double t2 = alpha_sx_ * (exp(-rc2*q2pG2) * q2pG2i );
              SumExpQpG2 += t1 + t2;
              if ( compute_stress )
              {
                const double tq1pGx = q1pG.x;
                const double tq1pGy = q1pG.y;
                const double tq1pGz = q1pG.z;

                const double tq2pGx = q2pG.x;
                const double tq2pGy = q2pG.y;
                const double tq2pGz = q2pG.z;

                // factor 2.0: derivative of (q+G)^2
                const double fac1 = t1 * 2.0 * ( rc2 + q1pG2i );
                const double fac2 = t2 * 2.0 * ( rc2 + q2pG2i );
                sigma_sumexp[0] += fac1 * tq1pGx * tq1pGx +
                                   fac2 * tq2pGx * tq2pGx;
                sigma_sumexp[1] += fac1 * tq1pGy * tq1pGy +
                                   fac2 * tq2pGy * tq2pGy;
                sigma_sumexp[2] += fac1 * tq1pGz * tq1pGz +
                                   fac2 * tq2pGz * tq2pGz;
                sigma_sumexp[3] += fac1 * tq1pGx * tq1pGy +
                                   fac2 * tq2pGx * tq2pGy;
                sigma_sumexp[4] += fac1 * tq1pGy * tq1pGz +
                                   fac2 * tq2pGy * tq2pGz;
                sigma_sumexp[5] += fac1 * tq1pGx * tq1pGz +
                                   fac2 * tq2pGx * tq2pGz;
              }
            }
          }

          // Add weighted contribution to numerical correction:
          // add the term sum_G exp(-a*|q+G|^2)/|q+G|^2 to the numerical
          // correction, if this is the first iKpoint.
          //
          // divide weight by 2 as we implicitly counted kpoint j and symmetric
          if ( iRotationStep==0 )
          {
            if ( iKpi==iKpj )
            {
              num_corr[iKpi] += SumExpQpG2 * 0.5 * wf.weight(iKpj);
              for ( int k = 0; k < 6; k++ )
              {
                sigma_num_corr[iKpi][k] += sigma_sumexp[k]*0.5*wf.weight(iKpi);
              }
            }
            else
            {
              num_corr[iKpi] += SumExpQpG2 * 0.5 * wf.weight(iKpj);
              num_corr[iKpj] += SumExpQpG2 * 0.5 * wf.weight(iKpi);
              for ( int k = 0; k < 6; k++ )
              {
                sigma_num_corr[iKpi][k] += sigma_sumexp[k]*0.5*wf.weight(iKpj);
                sigma_num_corr[iKpj][k] += sigma_sumexp[k]*0.5*wf.weight(iKpi);
              }
            }
          }

          // get occupation numbers for kpoint j
          const double* occ = sdj.occ_ptr();
          for ( int i = 0; i<sdj.nstloc(); i++ )
            occ_kj_[i]=occ[cj.jglobal(i)];

          // loop over the states at kpoint i
          for ( int i = 0; i < nStatesKpi_; i++ )
          {
            // loop over the states at kpoint j
            for ( int j = 0; j < sdj.nstloc(); j++ )
            {
              // check if something to compute for this pair
              if ( ( occ_ki_[i]!=0.0 && wf.weight(iKpi)!=0.0 )
                || ( occ_kj_[j]!=0.0 && wf.weight(iKpj)!=0.0 ) )
              {
                // compute the pair densities
                // contributions from both (ki,kj) and (ki,-kj)
                // rhor1: conj(psi(kj))*psi(ki)
                // rhor2: conj(psi(-kj))*psi(ki) == psi(kj)*psi(ki)
                for ( int ir = 0; ir < np012loc_; ir++ )
                {
                  rhor1_[ir] = conj( statej_[j][ir] ) * statei_[i][ir];
                  rhor2_[ir] =       statej_[j][ir]   * statei_[i][ir];
                }

                // Fourier transform the pair density to obtain rho(G).
                vft_->forward(&rhor1_[0], &rhog1_[0]);
                vft_->forward(&rhor2_[0], &rhog2_[0]);

                // initialize contrib of the states psi(kj,j) to the exchange
                // energy associated to state psi(ki,i)
                double ex_ki_i_kj_j = 0.0;

                double sigma_sum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

                for ( int ig = 0; ig < ngloc; ig++ )
                {
                  // Add the values of |rho1(G)|^2 * vint(q+G)
                  // to the exchange energy.
                  const double t1 = norm(rhog1_[ig]) * int_pot1_[ig];
                  const double t2 = norm(rhog2_[ig]) * int_pot2_[ig];
                  ex_ki_i_kj_j += t1;
                  ex_ki_i_kj_j += t2;

                  if ( compute_stress )
                  {
                    D3vector G(vbasis_->gx(ig+ngloc*0),
                               vbasis_->gx(ig+ngloc*1),
                               vbasis_->gx(ig+ngloc*2));
                    D3vector q1pG(q1+G);
                    D3vector q2pG(q2+G);
                    const double q1pG2 = norm2(q1pG);
                    const double q2pG2 = norm2(q2pG);
                    // dvint(g2) = d vint(g2)/d g2
                    const double d_int_pot1 = dvint(q1pG2);
                    const double d_int_pot2 = dvint(q2pG2);

                    const double tq1pGx = q1pG.x;
                    const double tq1pGy = q1pG.y;
                    const double tq1pGz = q1pG.z;

                    const double tq2pGx = q2pG.x;
                    const double tq2pGy = q2pG.y;
                    const double tq2pGz = q2pG.z;

                    // factor 2.0: derivative of (q+G)^2
                    const double fac1 = -2.0 * norm(rhog1_[ig]) * d_int_pot1;
                    const double fac2 = -2.0 * norm(rhog2_[ig]) * d_int_pot2;

                    sigma_sum[0] += fac1 * tq1pGx * tq1pGx;
                    sigma_sum[1] += fac1 * tq1pGy * tq1pGy;
                    sigma_sum[2] += fac1 * tq1pGz * tq1pGz;
                    sigma_sum[3] += fac1 * tq1pGx * tq1pGy;
                    sigma_sum[4] += fac1 * tq1pGy * tq1pGz;
                    sigma_sum[5] += fac1 * tq1pGx * tq1pGz;

                    sigma_sum[0] += fac2 * tq2pGx * tq2pGx;
                    sigma_sum[1] += fac2 * tq2pGy * tq2pGy;
                    sigma_sum[2] += fac2 * tq2pGz * tq2pGz;
                    sigma_sum[3] += fac2 * tq2pGx * tq2pGy;
                    sigma_sum[4] += fac2 * tq2pGy * tq2pGz;
                    sigma_sum[5] += fac2 * tq2pGx * tq2pGz;

                  }

                  if ( dwf )
                  {
                    // compute rhog1_[G]*V(|G+q1|) and rhog2_[G]*V(|G+q1|)
                    rhog1_[ig] *= int_pot1_[ig];
                    rhog2_[ig] *= int_pot2_[ig];
                  }
                }

                if ( dwf )
                {
                  // Backtransform rhog[G]*V(|q+G|)
                  vft_->backward(&rhog1_[0], &rhor1_[0]);
                  vft_->backward(&rhog2_[0], &rhor2_[0]);
                }

                // if iKpi=iKpj, add this contribution to the
                // exchange energy of state psi(ki,i)
                if ( iKpi==iKpj )
                {
                  // case iKpi=iKpj
                  // count only for psi(ki,i)
                  //
                  // compute the weights:
                  //
                  // => divide the weight of kpoint j by 2 as we implicitly
                  //    counted both kj and -kj in ex_ki_i_kj_j.
                  // => take into account the occupation of the state psi(kj,j)
                  const double exfac = - ( 4.0 * M_PI / omega ) * spinFactor;
                  double weight = exfac * 0.5 * wf.weight(iKpj) * occ_kj_[j];

                  // add contribution to exchange energy
                  double fac = weight * wf.weight(iKpi) * 0.5 * occ_ki_[i];
                  exchange_sum += fac * ex_ki_i_kj_j;

                  if ( compute_stress )
                  {
                    sigma_exhf_[0] += fac * (ex_ki_i_kj_j - sigma_sum[0]) /
                                      omega;
                    sigma_exhf_[1] += fac * (ex_ki_i_kj_j - sigma_sum[1]) /
                                      omega;
                    sigma_exhf_[2] += fac * (ex_ki_i_kj_j - sigma_sum[2]) /
                                      omega;
                    sigma_exhf_[3] += fac * ( -sigma_sum[3] ) / omega;
                    sigma_exhf_[4] += fac * ( -sigma_sum[4] ) / omega;
                    sigma_exhf_[5] += fac * ( -sigma_sum[5] ) / omega;
                  }

                  if (dwf)
                  {
                    // acumulate weighted contributions
                    // Psi_j,kj(r) * TF( rhog[G]/|q+G|^2 ) and symmetric
                    // in dpsi_i. We take now into account the mixing coeff
                    for ( int ir = 0; ir < np012loc_; ir++ )
                    {
                      dstatei_[i][ir] += ( statej_[j][ir] * rhor1_[ir] +
                        conj(statej_[j][ir] ) * rhor2_[ir]) * weight;
                    }
                  }
                }

                // if iKpi/=iKpj, add this contrib to the exch energy of both
                // states psi(ki,i) and psi(kj,j)
                // (as ex(ki,i,kj,j)=ex(kj,j,ki,i))
                // this way, we can avoid computing ex(ki,i,kj,j) for kj<ki.
                else
                {
                  // case iKpj>iKpi
                  // count for both psi(ki,i) and psi(kj,j)
                  //
                  // compute the weights:
                  //
                  // => divide the weight of kpoints j (resp. i) by 2 as we
                  //    implicitly counted both kj and -kj (resp ki and -ki)
                  //    in ex_ki_i_kj_j.
                  // => take into account the occupation of the state psi(kj,j)
                  //    (resp psi(ki,i))
                  double weighti = exfac * 0.5 * wf.weight(iKpi) * occ_ki_[i];
                  double weightj = exfac * 0.5 * wf.weight(iKpj) * occ_kj_[j];

                  // add contribution to exchange energy
                  double fac_i = weightj * wf.weight(iKpi) * 0.5 * occ_ki_[i];
                  double fac_j = weighti * wf.weight(iKpj) * 0.5 * occ_kj_[j];
                  double fac = fac_i + fac_j;
                  exchange_sum += fac * ex_ki_i_kj_j;

                  if ( compute_stress )
                  {
                    sigma_exhf_[0] += fac * (ex_ki_i_kj_j - sigma_sum[0]) /
                                      omega;
                    sigma_exhf_[1] += fac * (ex_ki_i_kj_j - sigma_sum[1]) /
                                      omega;
                    sigma_exhf_[2] += fac * (ex_ki_i_kj_j - sigma_sum[2]) /
                                      omega;
                    sigma_exhf_[3] += fac * ( -sigma_sum[3] ) / omega;
                    sigma_exhf_[4] += fac * ( -sigma_sum[4] ) / omega;
                    sigma_exhf_[5] += fac * ( -sigma_sum[5] ) / omega;
                  }

                  if (dwf)
                  {
                    // acumulate weighted contributions in dpsi_j and dpsi_j.
                    // the correspondances between rho12 and rho 21 are given by
                    //
                    //              /
                    // rho1_12(r) = | psi2*(r')*psi1(r')/|r-r'|dr'
                    //              /
                    //
                    //              /
                    // rho2_12(r) = | psi2(r')*psi1(r')/|r-r'|dr'
                    //              /
                    //
                    // => rho1_21  = rho1_12*
                    // => rho2_21  = rho2_12
                    //
                    // We take also into account the mixing coefficient.
                    //
                    for ( int ir = 0; ir < np012loc_; ir++ )
                    {
                      dstatei_[i][ir] += ( statej_[j][ir] * rhor1_[ir] +
                        conj( statej_[j][ir] ) * rhor2_[ir] ) * weightj;
                      dstatej_[j][ir] += ( statei_[i][ir] * conj(rhor1_[ir]) +
                        conj( statei_[i][ir] )* rhor2_[ir] ) * weighti;
                    }
                  }
                }
              } // something to do for pair (i,j)
            } // for j
          } // for i

          if (dwf)
          {
            // add the g space contributions to each state derivative of kpj
            for ( int i = 0; i < dsdj.nstloc(); i++ )
            {
              // compute the g space state derivative contribution
              wftj_->forward(&(dstatej_[i])[0],&buffer_dstate_[0]);

              // add the g the result to the state derivative of kpj
              complex<double> *p=dcj.valptr(i*dcj.mloc());
              for ( int j=0; j<dcj.mloc(); j++ )
                p[j]+=buffer_dstate_[j];
            }
          }
          delete wftj_;
        } // for iKpj

        // End of loop over kpoints j
        // finish to rotate the columns

        if (dwf)
        {
          CompleteReceivingForces(iRotationStep);

          CompleteSendingForces(iRotationStep);

          // add the g coordinate expression contributions to each
          // state derivative of kpi, store everything in the send buffer
          for ( int i=0; i<nStatesKpi_; i++ )
          {
            // transform contribution to g coordinates
            wfti_->forward(&(dstatei_[i])[0], &buffer_dstate_[0]);

            // sum up contributions in send buffer
            complex<double> *p1=&force_kpi_[i*dci.mloc()];
            complex<double> *p2=&send_buf_forces_[i*dci.mloc()];
            for ( int j=0; j<dci.mloc(); j++ )
              p2[j] = p1[j] + buffer_dstate_[j];
          }

          // Start forces permutation
          StartForcesPermutation(dci.mloc());
        }

        CompleteSendingOccupations(iRotationStep);

        // start occupations permutation
        StartOccupationsPermutation();

        // set the new number of local states
        nStatesKpi_ = nNextStatesKpi_;
      } // end iRotationStep

      //  end of rotation of the states of kpoint i from this point

      // wait for all communications to be completed
      {
        // complete all permutations except forces
        CompleteReceivingStates(1);
        CompleteSendingStates(1);
        CompleteReceivingOccupations(1);
        CompleteSendingOccupations(1);

        if (dwf)
        {
          // complete forces permutation
          CompleteReceivingForces(1);
          CompleteSendingForces(1);
        }

        // Terminate permutation
        FreePermutation();
      }

      // divergence corrections
      if ( alpha_sx_ != 0.0 )
      {
        const double integ = alpha_sx_ * 4.0 * M_PI * sqrt(M_PI) /
          ( 2.0 * rcut_ );
        const double vbz = pow(2.0*M_PI,3.0) / omega;

        for ( int i = 0; i < sdi.nstloc(); i++ )
        {
          double div_corr = 0.0;

          const double div_corr_1 = exfac * num_corr[iKpi] * occ_ki_[i];
          div_corr += div_corr_1;
          const double e_div_corr_1 = -0.5 * div_corr_1 * occ_ki_[i] *
                       wf.weight(iKpi);
          exchange_sum += e_div_corr_1;

          // contributions to stress from div_corr_1;
          const double fac = 0.5 * exfac * occ_ki_[i] * occ_ki_[i] *
                             wf.weight(iKpi);
          if ( compute_stress )
          {
            sigma_exhf_[0] += ( e_div_corr_1+ fac * sigma_num_corr[iKpi][0] ) /
                              omega;
            sigma_exhf_[1] += ( e_div_corr_1+ fac * sigma_num_corr[iKpi][1] ) /
                              omega;
            sigma_exhf_[2] += ( e_div_corr_1+ fac * sigma_num_corr[iKpi][2] ) /
                              omega;
            sigma_exhf_[3] += ( fac * sigma_num_corr[iKpi][3] ) / omega;
            sigma_exhf_[4] += ( fac * sigma_num_corr[iKpi][4] ) / omega;
            sigma_exhf_[5] += ( fac * sigma_num_corr[iKpi][5] ) / omega;
          }

          // rcut*rcut divergence correction
          if ( vbasis_->mype() == 0 )
          {
            const double div_corr_2 = - alpha_sx_ * exfac *
               rcut_ * rcut_ * occ_ki_[i] * wf.weight(iKpi);
            div_corr += div_corr_2;
            const double e_div_corr_2 = -0.5 * div_corr_2 * occ_ki_[i];
            exchange_sum += e_div_corr_2 * wf.weight(iKpi);

            // contributions of div_corr_2 to stress
            if ( compute_stress )
            {
              const double fac = wf.weight(iKpi);
              sigma_exhf_[0] += fac * e_div_corr_2 / omega;
              sigma_exhf_[1] += fac * e_div_corr_2 / omega;
              sigma_exhf_[2] += fac * e_div_corr_2 / omega;
            }

            const double div_corr_3 = - exfac * integ/vbz * occ_ki_[i];
            div_corr += div_corr_3;
            const double e_div_corr_3 = -0.5 * div_corr_3 * occ_ki_[i];
            exchange_sum += e_div_corr_3 * wf.weight(iKpi);
            // no contribution to stress from div_corr_3
          }

          // contribution of divergence corrections to forces on wave functions
          if (dwf)
          {
            // sum the partial contributions to the correction for state i
            gcontext_.dsum('C', 1, 1, &div_corr, 1);

            // add correction to the derivatives of state i
            complex<double> *ps=ci.valptr(i*ci.mloc());
            complex<double> *pf1=dci.valptr(i*dci.mloc());
            complex<double> *pf2=&force_kpi_[i*dci.mloc()];
            for ( int j = 0; j < dci.mloc(); j++ )
              pf1[j] += pf2[j] - ps[j] * div_corr;
          }
        } // for i
      } // if alpha_sx_
      delete wfti_;
    } // for iKpi
  } // for ispin

  // sum contributions to the exchange energy
  gcontext_.dsum(1, 1, &exchange_sum, 1);

  // sum stress tensor contributions
  if ( compute_stress )
    gcontext_.dsum(6,1,&sigma_exhf_[0],6);

  tm.stop();
#endif
  return exchange_sum;
}

////////////////////////////////////////////////////////////////////////////////
double ExchangeOperator::compute_exchange_at_gamma_(const Wavefunction &wf,
  Wavefunction* dwf, bool compute_stress)
{
  Timer tm;
  Timer tmb;

  wfc_ = wf;

  cout << setprecision(10);
  const double omega = wfc_.cell().volume();
  const int nspin = wfc_.nspin();

  // spin factor for the pair densities: 0.5 if 1 spin, 1 if nspin==2
  const double spinFactor=0.5*nspin;

  // total exchange energy
  double exchange_sum = 0.0;

  sigma_exhf_ = 0.0;
  const double *const g_x = vbasis_->gx_ptr(0);
  const double *const g_y = vbasis_->gx_ptr(1);
  const double *const g_z = vbasis_->gx_ptr(2);

  for ( int isp_loc = 0; isp_loc < wfc_.nsp_loc(); ++isp_loc )
  {
    const int ispin = wf.isp_global(isp_loc);

    SlaterDet& sd = *(wfc_.sd(isp_loc,0));
    ComplexMatrix& c = sd.c();
    const int nst = sd.nst();

    assert(wf.sd(isp_loc,0));
    // create Fourier transform object for wave functions
    if ( !wft_ )
    {
      if ( isp_loc >= 0 )
        wft_ = new FourierTransform(wf.sd(isp_loc,0)->basis(),
               np0v_,np1v_,np2v_);
    }

    // if using bisection, localize the wave functions
    if ( use_bisection_ )
    {
      tmb.reset();
      tmb.start();
      int maxsweep = 50;
      if ( s_.ctrl.debug.find("BISECTION_MAXSWEEP") != string::npos )
      {
        // override tolerance for bisection
        istringstream is(s_.ctrl.debug);
        string s;
        is >> s >> maxsweep;
        if ( MPIdata::onpe0() )
          cout << " override bisection maxsweep value: maxsweep = "
               << maxsweep << endl;
        assert(maxsweep >= 0);
      }
      double tol = 0.001;
      if ( s_.ctrl.debug.find("BISECTION_TOL") != string::npos )
      {
        // override tolerance for bisection
        istringstream is(s_.ctrl.debug);
        string s;
        is >> s >> tol;
        if ( MPIdata::onpe0() )
          cout << " override bisection tol value: tol = " << tol << endl;
        assert(tol > 0.0);
      }
#if TIMING
      Timer tmbtransf;
      tmbtransf.start();
#endif
      bisection_[isp_loc]->compute_transform(*wfc_.sd(isp_loc,0),maxsweep,tol);
#if TIMING
      tmbtransf.stop();
      Timer tmbcomploc;
      tmbcomploc.start();
#endif
      bisection_[isp_loc]->compute_localization(s_.ctrl.btHF);
#if TIMING
      tmbcomploc.stop();
#endif
      // copy of localization vector from Bisection object
      localization_ = bisection_[isp_loc]->localization();

      // copy the orthogonal transformation u to uc_[isp_loc]
      *uc_[isp_loc] = bisection_[isp_loc]->u();

      bool distribute = s_.ctrl.debug.find("BISECTION_NODIST") == string::npos;
      if ( distribute )
      {
        // define a permutation ordering states by increasing degree
        // permute states according to the order defined by the
        // localization vector

        // compute the degree of the vertices of the exchange graph
        // using the localization vector
#if TIMING
        Timer tmb_ov;
        tmb_ov.start();
#endif
        vector<int> degree(nst);
        for ( int i = 0; i < nst; i++ )
        {
          int count = 0;
          for ( int j = 0; j < nst; j++ )
          {
            if ( bisection_[isp_loc]->overlap(localization_,i,j) )
              count++;
          }
          degree[i] = count;
        }
#if TIMING
        tmb_ov.stop();
        if ( MPIdatat::onpe0() )
        {
          cout << setprecision(3);
          cout << " ExchangeOperator: bisection overlap time: "
             << tmb_ov.real() << " s" << endl;
        }
#endif

        // permutation index
        vector<int> index(nst);
        for ( int j = 0; j < index.size(); j++ )
          index[j] = j;

        // Create function object for comparison of degree
        VectorLess<int> degree_less(degree);
        sort(index.begin(), index.end(), degree_less);
        // At this point degree[index[i]] <= degree[index[j]] if i < j
        for ( int i = 0; i < index.size()-1; i++ )
          assert(degree[index[i]] <= degree[index[i+1]]);

#if DEBUG
        if ( MPIdata::onpe0() )
        {
          cout << "degree order after sort:" << endl;
          for ( int j = 0; j < index.size(); j++ )
            cout << j << " -> " << index[j]
                 << "  " << degree[index[j]]
                 << endl;
        }
#endif

        // distribute the states to process columns in round robin fashion
        // Assume that the states are initially ordered by increasing degree
        // i.e. degree(index_[i]) < degree(index_[j]) if i < j

        const int nb = uc_[isp_loc]->nb();

        vector<int> distrib_index(nst);
        int ibase = 0;
        int icol = 0;
        for ( int i = 0; i < nst; i++ )
        {
          // check if next slot is beyond n
          if ( ibase + icol * nb >= nst )
          {
            // restart in column 0 with incremented ibase
            icol = 0;
            ibase++;
          }
          distrib_index[ibase + icol * nb] = i;
          icol++;
        }

        // combine index[i] and distrib_index[i]
        vector<int> itmp(index.size());
        for ( int i = 0; i < index.size(); i++ )
          itmp[i] = index[distrib_index[i]];
        index = itmp;

#if DEBUG
        if ( MPIdata::onpe0() )
        {
          cout << "index after round robin distrib:" << endl;
          for ( int j = 0; j < index.size(); j++ )
             cout << j << " -> " << index[j] << endl;
        }
#endif
        // apply the permutation defined by index to the localization vector
        vector<long int> loc_tmp(localization_.size());
        for ( int i = 0; i < index.size(); i++ )
          loc_tmp[i] = localization_[index[i]];
        localization_ = loc_tmp;

        // apply the permutation defined by index to the occupation vector
        vector<double> occ_tmp(nst);
        for ( int i = 0; i < index.size(); i++ )
          occ_tmp[i] = sd.occ(index[i]);
        sd.set_occ(occ_tmp);

        // compute a pivot vector containing a sequence of transpositions
        // equivalent to the permutation defined by index[i]
        // compute the inverse of index
        vector<int> index_inv(index.size());
        for ( int i = 0; i < index.size(); i++ )
          index_inv[index[i]] = i;
        assert(index_inv.size()==index.size());

        vector<int> pivot;
        for ( int i = 0; i < index.size(); i++ )
        {
          int j = index_inv[i];

          int tmp = index[i];
          index[i] = index[j];
          index[j] = tmp;

          // update elements of index_inv
          index_inv[index[i]] = i;
          index_inv[index[j]] = j;

          pivot.push_back(j);
        }
        assert(pivot.size()==index.size());
#if DEBUG
        if ( MPIdata::onpe0() )
        {
          cout << "pivot:" << endl;
          for ( int j = 0; j < pivot.size(); j++ )
            cout << j << " -> " << pivot[j] << endl;
        }
#endif

        // create a local pivot vector on this process (size uc->nloc())
        // this vector must be replicated on all tasks of the
        // process grid columns
        const int nloc = uc_[isp_loc]->nloc();
        vector<int> locpivot(nloc);
        // fill the local pivot vector on all tasks
        // add 1 to index values for lapack fortran index convention
        for ( int j=0; j < nloc; j++ )
        {
          int jglobal = uc_[isp_loc]->jglobal(j);
          locpivot[j] = pivot[jglobal]+1;
        }
#if 0
        for ( int ipe = 0; ipe < uc.context().size(); ipe++ )
        {
          uc.context().barrier();
          if ( ipe == uc.context().mype() )
          {
            cout << "locpivot:" << endl;
            for ( int j = 0; j < locpivot.size(); j++ )
              cout << ipe << ": " << j << " -> " << locpivot[j] << endl;
          }
        }
#endif

#if 0
        if ( u_->context().size() == 1 )
        {
          // cout << " uc before perm: " << endl;
          // cout << uc;
          // local permutation
          assert(locpivot.size()==u_->n());
          double *p = uc.valptr(0);
          const int mloc = uc.mloc();
          for ( int i = locpivot.size()-1; i >= 0; i-- )
          {
            const int j = locpivot[i]-1;
            // swap columns i and j of u_->c()
            cout << " swap " << i << " " << j << endl;
            for ( int k = 0; k < mloc; k++ )
            {
              double tmp = p[i*mloc+k];
              p[i*mloc+k] = p[j*mloc+k];
              p[j*mloc+k] = tmp;
            }
          }
          // cout << " uc after perm: " << endl;
          // cout << uc;
        }
#endif

#if TIMING
        Timer tmblapiv;
        tmblapiv.start();
#endif
        // apply the permutation to the columns of uc
        uc_[isp_loc]->lapiv('B','C',&locpivot[0]);
#if TIMING
        tmblapiv.stop();

        if ( MPIdata::sd_rank() == 0 )
        {
          cout << setprecision(3);
          cout << " ExchangeOperator: bisection size time: "
             << tmbsize.real() << " s" << endl;
          cout << setprecision(3);
          cout << " ExchangeOperator: bisection pair time: "
             << tmbpair.real() << " s" << endl;
          cout << setprecision(3);
          cout << " ExchangeOperator: bisection lapiv time: "
             << tmblapiv.real() << " s" << endl;
        }
#endif

#if DEBUG
        // recompute the degree of the vertices of the exchange graph
        for ( int i = 0; i < nst; i++ )
        {
          int count = 0;
          for ( int j = 0; j < nst; j++ )
          {
            if ( bisection_[isp_loc]->overlap(localization_,i,j) )
              count++;
          }
          degree[i] = count;
        }

        if ( MPIdata::sd_rank() == 0 )
        {
          cout << "degree after permutation:" << endl;
          for ( int j = 0; j < degree.size(); j++ )
            cout << " deg[" << j << "] = " << degree[j] << endl;
        }

        // print localization vectors and overlaps after distribution
        if ( MPIdata::sd_rank() == 0 )
        {
          int sum = 0;
          for ( int i = 0; i < nst; i++ )
          {
            int count = 0;
            for ( int j = 0; j < nst; j++ )
            {
              if ( bisection_[isp_loc]->overlap(localization_,i,j) )
                count++;
            }
            cout << "localization[" << i << "]: "
                 << localization_[i] << " "
                 << bitset<30>(localization_[i])
                 << "  overlaps: "
                 << count << endl;
            sum += count;
          }
          cout << " Bisection::localize: total overlaps: " << sum << " / "
               << nst*nst << " = " << ((double) sum)/(nst*nst) << endl;
        }
#endif

      }
      else
      {
        if ( MPIdata::onpe0() )
          cout << " ExchangeOperator: bisection distribution disabled" << endl;
      } // if distribute

#if TIMING
      Timer tmbfwd;
      tmbfwd.start();
#endif
      bisection_[isp_loc]->forward(*uc_[isp_loc], *wfc_.sd(isp_loc,0));
#if TIMING
      tmbfwd.stop();
#endif

      tmb.stop();
    } // if use_bisection_

    tm.start();

    // compute exchange

    // real-space local states -> statej_[i][ir]
    // loop over states 2 by 2
    for ( int i = 0; i < sd.nstloc()-1; i+=2 )
    {
      wft_->backward(c.cvalptr(i*c.mloc()), c.cvalptr((i+1)*c.mloc()),&tmp_[0]);
      double *p = (double *)&tmp_[0];
#pragma omp parallel for
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        statej_[i][ir]=p[2*ir];
        statej_[i+1][ir]=p[2*ir+1];
      }
    }

    if ( sd.nstloc() % 2 == 1 )
    {
      int i = sd.nstloc()-1;
      wft_->backward(c.cvalptr(i*c.mloc()), &(statej_[i])[0]);
    }

    SlaterDet& dsd = *(dwf->sd(isp_loc,0));
    ComplexMatrix& dc = dsd.c();

    if (dwf)
    {
      // reset real space derivatives dstatej_[i][ir]
      for ( int i = 0; i < dsd.nstloc(); i++ )
        for ( int j = 0; j < np012loc_; j++ )
          dstatej_[i][j] = 0.0;
    }

    const int ngloc = vbasis_->localsize();
    // correction term: sum_(G) exp(-rcut_^2*|G|^2)/|G|^2
    const double exfac = - ( 4.0 * M_PI / omega ) * spinFactor;
    double SumExpG2 = 0.0;
    double sigma_sumexp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const double *g2 = vbasis_->g2_ptr();
    const double *g2i = vbasis_->g2i_ptr();
    const double rc2 = rcut_*rcut_;

    // divergence correction (Coulomb potential part only)
    // The interaction potential is
    // beta_sx*(1/r) + (alpha_sx-beta_sx)*erf(mu*r)/r
    // The coefficient of the long range term is alpha_sx
    // subtract alpha_sx * exp(-rc2*G^2)/G^2
    if ( alpha_sx_ != 0.0 )
    {
      for ( int ig = 0; ig < ngloc; ig++ )
      {
        // factor 2.0: real basis
        const double tg2i = g2i[ig];
        double t = alpha_sx_ * 2.0 * exp( - rc2 * g2[ig] ) * tg2i;
        SumExpG2 += t;

        if ( compute_stress )
        {
          const double tgx = g_x[ig];
          const double tgy = g_y[ig];
          const double tgz = g_z[ig];
          // factor 2.0: derivative of G^2
          const double fac = t * 2.0 * ( rc2 + tg2i );
          sigma_sumexp[0] += fac * tgx * tgx;
          sigma_sumexp[1] += fac * tgy * tgy;
          sigma_sumexp[2] += fac * tgz * tgz;
          sigma_sumexp[3] += fac * tgx * tgy;
          sigma_sumexp[4] += fac * tgy * tgz;
          sigma_sumexp[5] += fac * tgz * tgx;
        }
      }
    }

    // local occupation numbers
    const double* occ = sd.occ_ptr();
    for ( int i = 0; i < sd.nstloc(); i++ )
      occ_kj_[i]=occ[c.jglobal(i)];

    // number of states to be sent
    nStatesKpi_ = sd.nstloc();

    // occupation numbers of circulating states
    for ( int i = 0; i < nStatesKpi_; i++ )
      occ_ki_[i] = occ_kj_[i];

    // copy local states into circulating states
    const complex<double> *p = c.cvalptr(0);
    for ( int i = 0; i < nStatesKpi_ * c.mloc(); i++ )
      state_kpi_[i]=p[i];

    // initialize circulating derivatives
    if (dwf)
    {
      for ( int i = 0; i < nStatesKpi_ * dc.mloc(); i++ )
        force_kpi_[i]=0.0;
    }

    // initiate send nStatesKpi_ and receive nNextStatesKpi_
    InitPermutation();

#if LOAD_MATRIX
    // collect number of processed pairs in array load_matrix
    vector<int> load_matrix(MPIdata::sd_size(),0);
#endif

    // Start rotation of circulating states
    for ( int iRotationStep = 0; iRotationStep < MPIdata::nstb();
          iRotationStep++ )
    {
      // generate a list of pairs of overlapping states
      int nPair = 0;
      vector<int> first_member_of_pair;
      vector<int> second_member_of_pair;

      // flag indicating that a circulating state has been used
      vector<int> useState(nStatesKpi_,0);

      // finish receiving occupations in occ_ki_[]
      CompleteReceivingOccupations(iRotationStep);

      // loop over circulating states
      for ( int i = 0; i < nStatesKpi_; i++ )
      {
        // original column index of circulating state i
        int iColI = MPIdata::istb() - iRotationStep;
        iColI = ( iColI < 0 ) ? iColI + MPIdata::nstb() : iColI;

        // global index of circulating state i
        int iGlobI = c.jglobal(iColI,i);

        // loop over fixed states
        for ( int j = 0; j < sd.nstloc(); j++ )
        {
          // check if there is something to compute for this pair
          if ( occ_ki_[i]!=0.0 || occ_kj_[j]!=0.0 )
          {
            // global index of fixed state j
            int iGlobJ = c.jglobal(j);

            // determine the overlap between those two states
            bool overlap_ij = ( !use_bisection_ ||
              bisection_[isp_loc]->overlap(localization_,iGlobI,iGlobJ) );

            // use the chess board condition to
            // optimize the distribution of work on
            // each process:
            // - if i and j have different parity, use
            // the condition i<j
            // - if i and j have same parity, use
            // the condition i>=j
            int parity_i = iGlobI & 1;
            int parity_j = iGlobJ & 1;

            if ( parity_i == parity_j )
            {
              if ( iGlobI >= iGlobJ && overlap_ij )
              {
                first_member_of_pair.push_back( i );
                second_member_of_pair.push_back( j );
                nPair++;

                // circulating state i is used
                useState[i] = 1;
              }
            }
            else
            {
              if ( iGlobI < iGlobJ && overlap_ij )
              {
                first_member_of_pair.push_back( i );
                second_member_of_pair.push_back( j );
                nPair++;

                // circulating state i is used
                useState[i] = 1;
              }
            }
          }
        }
      }
      // pair list is complete

#if LOAD_MATRIX
      // collect nPair statistics if on row 0
      if ( MPIdata::g_rank() == 0 )
        load_matrix[iRotationStep*MPIdata::sd_size()+MPIdata::istb()] =
          nPair;
#endif

      // complete receiving states
      // note: this does nothing if iRotationStep == 0
      CompleteReceivingStates(iRotationStep);
      // circulating states in state_kpi_[i+j*mloc] can now be used

      // compute real space circulating states
      if ( nPair > 0 )
      {
        int i = 0;
        int j = 1;

        // try to associate states 2 by 2
        while ( i < nStatesKpi_ )
        {
          // seek the first next used state
          while ( i < nStatesKpi_ && !useState[i] ) i++;

          // seek the second next used state
          j=i+1;
          while ( j < nStatesKpi_ && !useState[j] ) j++;

          // if i is a valid index
          if ( i < nStatesKpi_ )
          {
            // if j is a valid index
            if ( j < nStatesKpi_ )
            {
              // back transform the couple (i,j)
              wft_->backward(&state_kpi_[i*c.mloc()],
                &state_kpi_[j*c.mloc()], &tmp_[0]);

              // copy the result in state[i] and state[j]
              double *p = (double *)&tmp_[0];
#pragma omp parallel for
              for ( int ir = 0; ir < np012loc_; ir++ )
              {
                statei_[i][ir]=p[2*ir];
                statei_[j][ir]=p[2*ir+1];
              }
            }

            // else, if there is only one state to transform
            else
            {
              wft_->backward(&state_kpi_[i*c.mloc()], &(statei_[i])[0]);
            }
          }

          // increment indices
          i = j + 1;
        }
      }

      // finish sending states in send_buf_states_
      CompleteSendingStates(iRotationStep);
      // send_buf_states_ can now be reused

      // copy the states to be sent in the send buffer
      for ( int i = 0; i < nStatesKpi_ * c.mloc(); i++ )
        send_buf_states_[i] = state_kpi_[i];

      if (dwf)
      {
        for ( int i = 0; i < nStatesKpi_; i++ )
          for ( int j = 0; j < np012loc_; j++ )
            dstatei_[i][j] = 0.0;
      }

      // nNextStatesKpi: number of states of next permutation step
      SetNextPermutationStateNumber();
      // start sending states in send_buf_states_
      StartStatesPermutation(c.mloc());

      // loop over pairs 2 by 2
      if ( nPair > 0 )
      {
        double ex_sum_1, ex_sum_2;
        double sigma_sum_1[6], sigma_sum_2[6];
        int iPair;
        for ( iPair=0; iPair<first_member_of_pair.size()-1; iPair+=2 )
        {
          int i1 = first_member_of_pair[iPair];
          int j1 = second_member_of_pair[iPair];
          int i2 = first_member_of_pair[iPair+1];
          int j2 = second_member_of_pair[iPair+1];

          // compute the pair densities
          // rhor = conjg(statei_(r)) * statej_(r)
          // note: gamma point, densities are real
          {
            // rhor1_ = psi_i1 * psi_j1 + i * psi_i2 * psi_j2
            double *p   = (double *)&rhor1_[0];
            double *pi1 = (double *)&statei_[i1][0];
            double *pi2 = (double *)&statei_[i2][0];
            double *pj1 = (double *)&statej_[j1][0];
            double *pj2 = (double *)&statej_[j2][0];

#pragma omp parallel for
            for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
            {
              p[ip]   = pi1[ip] * pj1[ip];
              p[ip+1] = pi2[ip] * pj2[ip];
            }
          }

          // Fourier transform the pair density
          vft_->forward(&rhor1_[0], &rhog1_[0], &rhog2_[0]);

          // compute contributions to the exchange energy and forces on wfs
          ex_sum_1 = 0.0;
          ex_sum_2 = 0.0;
          if ( compute_stress )
          {
            for ( int i = 0; i < 6; i++ )
            {
              sigma_sum_1[i] = 0.0;
              sigma_sum_2[i] = 0.0;
            }
          }

          for ( int ig = 0; ig < ngloc; ig++ )
          {
            // Add the values of |rho1(G)|^2*V(|G+q1|)
            // and |rho2(G)|^2*V(|G+q2|) to the exchange energy.
            // factor 2.0: real basis
            const double int_pot = vint(g2[ig]);
            const double t1 = 2.0 * norm(rhog1_[ig]) * int_pot;
            const double t2 = 2.0 * norm(rhog2_[ig]) * int_pot;
            ex_sum_1 += t1;
            ex_sum_2 += t2;

            if ( compute_stress )
            {
              // dvint(g2) = d vint(g2)/d g2
              const double d_int_pot = dvint(g2[ig]);
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              // factor 4.0: derivative of G^2 and real basis
              const double fac1 = -4.0 * norm(rhog1_[ig]) * d_int_pot;
              sigma_sum_1[0] += fac1 * tgx * tgx;
              sigma_sum_1[1] += fac1 * tgy * tgy;
              sigma_sum_1[2] += fac1 * tgz * tgz;
              sigma_sum_1[3] += fac1 * tgx * tgy;
              sigma_sum_1[4] += fac1 * tgy * tgz;
              sigma_sum_1[5] += fac1 * tgz * tgx;

              const double fac2 = -4.0 * norm(rhog2_[ig]) * d_int_pot;
              sigma_sum_2[0] += fac2 * tgx * tgx;
              sigma_sum_2[1] += fac2 * tgy * tgy;
              sigma_sum_2[2] += fac2 * tgz * tgz;
              sigma_sum_2[3] += fac2 * tgx * tgy;
              sigma_sum_2[4] += fac2 * tgy * tgz;
              sigma_sum_2[5] += fac2 * tgz * tgx;
            }

            if (dwf)
            {
              // compute rhog1_[G]*V(G) and rhog2_[G]*V(G)
              rhog1_[ig] *= int_pot;
              rhog2_[ig] *= int_pot;
            }
          }

          if (dwf)
          {
            // Backtransform rhog[G]/|q+G|^2
            vft_->backward(&rhog1_[0], &rhog2_[0], &rhor1_[0]);
          }

          // accumulate contributions to the exchange energy
          // first pair: (i1,j1)
          const double fac1 = 0.5 * exfac * occ_ki_[i1] * occ_kj_[j1];
          if ( ( i1==j1 ) && ( iRotationStep==0 ) )
          {
            exchange_sum += fac1 * ex_sum_1;

            if ( compute_stress )
            {
              sigma_exhf_[0] += fac1 * (ex_sum_1 - sigma_sum_1[0]) / omega;
              sigma_exhf_[1] += fac1 * (ex_sum_1 - sigma_sum_1[1]) / omega;
              sigma_exhf_[2] += fac1 * (ex_sum_1 - sigma_sum_1[2]) / omega;
              sigma_exhf_[3] += fac1 * ( -sigma_sum_1[3] ) / omega;
              sigma_exhf_[4] += fac1 * ( -sigma_sum_1[4] ) / omega;
              sigma_exhf_[5] += fac1 * ( -sigma_sum_1[5] ) / omega;
            }

            if (dwf)
            {
              const double weight = exfac * occ_kj_[j1];
              double *dp = (double *) &dstatei_[i1][0];
              double *pj = (double *) &statej_[j1][0];
              double *pr = (double *) &rhor1_[0];
#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
                dp[ip] += pj[ip] * pr[ip] * weight;
            }
          }
          else
          {
            exchange_sum += 2.0 * fac1 * ex_sum_1;

            sigma_exhf_[0] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[0]) / omega;
            sigma_exhf_[1] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[1]) / omega;
            sigma_exhf_[2] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[2]) / omega;
            sigma_exhf_[3] += 2.0 * fac1 * ( -sigma_sum_1[3] ) / omega;
            sigma_exhf_[4] += 2.0 * fac1 * ( -sigma_sum_1[4] ) / omega;
            sigma_exhf_[5] += 2.0 * fac1 * ( -sigma_sum_1[5] ) / omega;

            if (dwf)
            {
              double weighti = exfac * occ_ki_[i1];
              double weightj = exfac * occ_kj_[j1];
              double *dpi = (double *) &dstatei_[i1][0];
              double *dpj = (double *) &dstatej_[j1][0];
              double *pi = (double *) &statei_[i1][0];
              double *pj = (double *) &statej_[j1][0];
              double *pr = (double *) &rhor1_[0];

#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
              {
                dpi[ip] += pj[ip] * pr[ip] * weightj;
                dpj[ip] += pi[ip] * pr[ip] * weighti;
              }
            }
          }

          // second pair
          const double fac2 = 0.5 * exfac * occ_ki_[i2] * occ_kj_[j2];
          if ( ( i2==j2 ) && ( iRotationStep==0 ) )
          {
            exchange_sum += fac2 * ex_sum_2;

            sigma_exhf_[0] += fac2 * (ex_sum_2 - sigma_sum_2[0]) / omega;
            sigma_exhf_[1] += fac2 * (ex_sum_2 - sigma_sum_2[1]) / omega;
            sigma_exhf_[2] += fac2 * (ex_sum_2 - sigma_sum_2[2]) / omega;
            sigma_exhf_[3] += fac2 * ( -sigma_sum_2[3] ) / omega;
            sigma_exhf_[4] += fac2 * ( -sigma_sum_2[4] ) / omega;
            sigma_exhf_[5] += fac2 * ( -sigma_sum_2[5] ) / omega;

            if (dwf)
            {
              const double weight = exfac * occ_kj_[j2];
              double *dp = (double *) &dstatei_[i2][0];
              double *pj = (double *) &statej_[j2][0];
              double *pr = (double *) &rhor1_[0];
              pr = pr + 1;

#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
                dp[ip] += pj[ip] * pr[ip] * weight;
            }
          }
          else
          {
            exchange_sum += 2.0 * fac2 * ex_sum_2;

            sigma_exhf_[0] += 2.0 * fac2 * (ex_sum_2 - sigma_sum_2[0]) / omega;
            sigma_exhf_[1] += 2.0 * fac2 * (ex_sum_2 - sigma_sum_2[1]) / omega;
            sigma_exhf_[2] += 2.0 * fac2 * (ex_sum_2 - sigma_sum_2[2]) / omega;
            sigma_exhf_[3] += 2.0 * fac2 * ( -sigma_sum_2[3] ) / omega;
            sigma_exhf_[4] += 2.0 * fac2 * ( -sigma_sum_2[4] ) / omega;
            sigma_exhf_[5] += 2.0 * fac2 * ( -sigma_sum_2[5] ) / omega;

            if (dwf)
            {
              double weighti = exfac * occ_ki_[i2];
              double weightj = exfac * occ_kj_[j2];
              double *dpi = (double *) &dstatei_[i2][0];
              double *dpj = (double *) &dstatej_[j2][0];
              double *pi = (double *) &statei_[i2][0];
              double *pj = (double *) &statej_[j2][0];
              double *pr = (double *) &rhor1_[0];
              pr = pr + 1;

#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
              {
                dpi[ip] += pj[ip] * pr[ip] * weightj;
                dpj[ip] += pi[ip] * pr[ip] * weighti;
              }
            }
          }
        } // iPair

        // last pair if necessary
        if ( iPair < first_member_of_pair.size() )
        {
          int i1 = first_member_of_pair[iPair];
          int j1 = second_member_of_pair[iPair];

          // rhor1_ = psi_i1 * psi_j1
          double *p   = (double *)&rhor1_[0];
          double *pi1 = (double *)&statei_[i1][0];
          double *pj1 = (double *)&statej_[j1][0];
#pragma omp parallel for
          for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
          {
            p[ip]   = pi1[ip] * pj1[ip];
            p[ip+1] = 0.0;
          }

          vft_->forward(&rhor1_[0], &rhog1_[0]);

          ex_sum_1 = 0.0;
          if ( compute_stress )
          {
            for ( int i = 0; i < 6; i++ )
              sigma_sum_1[i] = 0.0;
          }
          for ( int ig = 0; ig < ngloc; ig++ )
          {
            // Add the values of |rho1(G)|^2*V(|G|)
            // and |rho2(G)|^2*V(|G|) to the exchange energy.
            // factor 2.0: real basis
            const double int_pot = vint(g2[ig]);
            const double t1 = 2.0 * norm(rhog1_[ig]) * int_pot;
            ex_sum_1 += t1;

            if (dwf)
            {
              rhog1_[ig] *= int_pot;
            }

            if ( compute_stress )
            {
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              // factor 4.0: derivative of G^2 and real basis
              const double fac = -4.0 * norm(rhog1_[ig]) * dvint(g2[ig]);
              sigma_sum_1[0] += fac * tgx * tgx;
              sigma_sum_1[1] += fac * tgy * tgy;
              sigma_sum_1[2] += fac * tgz * tgz;
              sigma_sum_1[3] += fac * tgx * tgy;
              sigma_sum_1[4] += fac * tgy * tgz;
              sigma_sum_1[5] += fac * tgz * tgx;
            }
          }

          if (dwf)
          {
            // Backtransform rhog[G]/|q+G|^2
            vft_->backward(&rhog1_[0], &rhor1_[0]);
          }

          const double fac1 = 0.5 * exfac * occ_ki_[i1] * occ_kj_[j1];
          if ( ( i1==j1 ) && ( iRotationStep==0 ) )
          {
            exchange_sum += fac1 * ex_sum_1;

            sigma_exhf_[0] += fac1 * (ex_sum_1 - sigma_sum_1[0]) / omega;
            sigma_exhf_[1] += fac1 * (ex_sum_1 - sigma_sum_1[1]) / omega;
            sigma_exhf_[2] += fac1 * (ex_sum_1 - sigma_sum_1[2]) / omega;
            sigma_exhf_[3] += fac1 * ( -sigma_sum_1[3] ) / omega;
            sigma_exhf_[4] += fac1 * ( -sigma_sum_1[4] ) / omega;
            sigma_exhf_[5] += fac1 * ( -sigma_sum_1[5] ) / omega;

            if (dwf)
            {
              const double weight = exfac * occ_kj_[j1];
              double *dp = (double *) &dstatei_[i1][0];
              double *pj = (double *) &statej_[j1][0];
              double *pr = (double *) &rhor1_[0];
#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
                dp[ip] += pj[ip] * pr[ip] * weight;
            }
          }
          else
          {
            exchange_sum += 2.0 * fac1 * ex_sum_1;

            sigma_exhf_[0] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[0]) / omega;
            sigma_exhf_[1] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[1]) / omega;
            sigma_exhf_[2] += 2.0 * fac1 * (ex_sum_1 - sigma_sum_1[2]) / omega;
            sigma_exhf_[3] += 2.0 * fac1 * ( -sigma_sum_1[3] ) / omega;
            sigma_exhf_[4] += 2.0 * fac1 * ( -sigma_sum_1[4] ) / omega;
            sigma_exhf_[5] += 2.0 * fac1 * ( -sigma_sum_1[5] ) / omega;

            if (dwf)
            {
              double weighti = exfac * occ_ki_[i1];
              double weightj = exfac * occ_kj_[j1];
              double *dpi = (double *) &dstatei_[i1][0];
              double *dpj = (double *) &dstatej_[j1][0];
              double *pi = (double *) &statei_[i1][0];
              double *pj = (double *) &statej_[j1][0];
              double *pr = (double *) &rhor1_[0];

#pragma omp parallel for
              for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
              {
                dpi[ip] += pj[ip] * pr[ip] * weightj;
                dpj[ip] += pi[ip] * pr[ip] * weighti;
              }
            }
          }
        } // if last pair
      } // if nPair > 0
      // End of loop over pairs

      if (dwf)
      {
        // finish receiving forces in force_kpi_[]
        CompleteReceivingForces(iRotationStep);
        // finish sending forces in send_buf_forces_[]
        CompleteSendingForces(iRotationStep);

        // add locally computed contributions to circulated forces
        {
          // copy force_kpi to to send_buf_forces
          for ( int i = 0; i < nStatesKpi_; i++ )
          {
            complex<double> *ps = &send_buf_forces_[i*dc.mloc()];
            complex<double> *pf = &force_kpi_[i*dc.mloc()];
            for ( int j = 0; j < dc.mloc(); j++)
              ps[j] = pf[j];
          }
          // backtransform computed forces to G coordinates and
          // add to send buffer
          if ( nPair > 0 )
          {
            int i = 0;
            int j = 1;

            // try to associate states 2 by 2
            while ( i < nStatesKpi_ )
            {
              // find the first next used state
              while ( i < nStatesKpi_ && !useState[i] ) i++;

              // find the second next used state
              j = i+1;
              while ( j < nStatesKpi_ && !useState[j] ) j++;

              // if i is a valid index
              if ( i < nStatesKpi_ )
              {
                // if j is a valid index
                if ( j < nStatesKpi_ )
                {
                  // store dstate[i] + i * dstate[j] in tmp_
                  double *p = (double *)&tmp_[0];
                  double *dpr = (double *)&dstatei_[i][0];
                  double *dpi = (double *)&dstatei_[j][0];
#pragma omp parallel for
                  for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
                  {
                    p[ip]=dpr[ip];
                    p[ip+1]=dpi[ip];
                  }

                  // transform the pair of states
                  wft_->forward( &(tmp_)[0], &buffer_forces_1_[0],
                    &buffer_forces_2_[0] );

                  // accumulate contributions in send buffer
                  complex<double> *ps1=&send_buf_forces_[i*dc.mloc()];
                  complex<double> *ps2=&send_buf_forces_[j*dc.mloc()];
                  for ( int k = 0; k < dc.mloc(); k++)
                  {
                    ps1[k] += buffer_forces_1_[k];
                    ps2[k] += buffer_forces_2_[k];
                  }
                }
                else
                {
                  // there is only one force to transform
                  wft_->forward(&(dstatei_[i])[0], &buffer_forces_1_[0]);

                  // accumulate contributions in send buffer
                  complex<double> *ps1=&send_buf_forces_[i*dc.mloc()];
                  for ( int k = 0; k < dc.mloc(); k++ )
                    ps1[k] += buffer_forces_1_[k];
                }
                i = j+1;
              }
            }
          }
        }
        // end back transform and addition of locally computed
        // forces to the send buffer
        StartForcesPermutation(dc.mloc());
      } // if dwf

      CompleteSendingOccupations(iRotationStep);
      StartOccupationsPermutation();

      // set the new number of local states
      nStatesKpi_ = nNextStatesKpi_;
    } // iRotationStep
    // end of rotation of the states of kpoint i from this point

#if LOAD_MATRIX
    // collect load_matrix
    vector<int> tmp_matrix(load_matrix);
    MPI_Allreduce(&tmp_matrix,&load_matrix,load_matrix.size(),MPI_INT,
      MPIdata::st_comm());
    if ( MPIdata::sd_rank() == 0 )
    {
      cout << " ExchangeOperator: load_matrix" << endl;
      const int nst = s_.wfc_.nst(isp_loc);
      int spreadsum = 0;
      for ( int irot = 0; irot < MPIdata::nstb(); irot++ )
      {
        int rowsum = 0;
        int wmin = nst*nst;
        int wmax = 0;
        for ( int icol = 0; icol < MPIdata::nstb(); icol++ )
        {
          int w = load_matrix[irot*MPIdata::nstb()+icol];
          cout << " " << setw(5) << w;
          rowsum += w;
          wmin = min(wmin,w);
          wmax = max(wmax,w);
        }
        int spread = abs(wmin-wmax);
        cout << "   " << setw(5) << rowsum
             << "   spread: " << spread << endl;
        spreadsum += spread;
      }
      cout << endl;
      // print sums of columns
      int rowcolsum = 0;
      for ( int icol = 0; icol < MPIdata::nstb(); icol++ )
      {
        int colsum = 0;
        for ( int irot = 0; irot < MPIdata::nstb(); irot++ )
        {
          int w = load_matrix[irot*MPIdata::nstb()+icol];
          colsum += w;
        }
        cout << " " << setw(5) << colsum;
        rowcolsum += colsum;
      }
      cout << "   " << setw(5) << rowcolsum
           << "   spread: " << spreadsum << endl;
      cout << " pair fraction: " << ((double) rowcolsum)/(0.5*nst*(nst+1))
           << endl;
    }
#endif

    // wait for all communications to be completed
    // complete all permutations except forces
    CompleteReceivingStates(1);
    CompleteSendingStates(1);
    CompleteReceivingOccupations(1);
    CompleteSendingOccupations(1);

    if (dwf)
    {
      // complete forces permutation
      CompleteReceivingForces(1);
      CompleteSendingForces(1);
    }

    FreePermutation();

    // transform accumulated real-space forces to G space
    // loop over pairs of states
    for ( int i = 0; i < nStatesKpi_-1; i+=2 )
    {
      double *p = (double *)&tmp_[0];
      double *dpr = (double *)&dstatej_[i][0];
      double *dpi = (double *)&dstatej_[i+1][0];
#pragma omp parallel for
      for ( int ip = 0; ip < 2*np012loc_; ip+=2 )
      {
        p[ip]=dpr[ip];
        p[ip+1]=dpi[ip];
      }
      // transform the pair of forces
      wft_->forward(&(tmp_)[0], &buffer_forces_1_[0], &buffer_forces_2_[0]);

      // accumulate contributions into dc
      complex<double> *p1=dc.valptr(i*dc.mloc());
      complex<double> *p2=dc.valptr((i+1)*dc.mloc());
      complex<double> *pf1=&force_kpi_[i*dc.mloc()];
      complex<double> *pf2=&force_kpi_[(i+1)*dc.mloc()];
      for ( int j = 0; j < dc.mloc(); j++ )
      {
        p1[j] = buffer_forces_1_[j] + pf1[j];
        p2[j] = buffer_forces_2_[j] + pf2[j];
      }
    }

    // transform the last state if necessary
    if ( nStatesKpi_ % 2 == 1 )
    {
      int i = nStatesKpi_ - 1;
      // transform the force
      wft_->forward(&(dstatej_[i])[0], &buffer_forces_1_[0]);

      // accumulate contributions into dc
      complex<double> *p1=dc.valptr(i*dc.mloc());
      complex<double> *pf1=&force_kpi_[i*dc.mloc()];
      for ( int j = 0; j < dc.mloc(); j++ )
        p1[j] = buffer_forces_1_[j] + pf1[j];
    }
    // dc now contains the forces

    // divergence corrections from long range Coulomb part
    if ( alpha_sx_ != 0.0 )
    {
      const double integ = alpha_sx_ * 4.0 * M_PI * sqrt(M_PI) /
        ( 2.0 * rcut_ );
      const double vbz = pow(2.0*M_PI,3.0) / omega;

      // correct the energy of state i
      for ( int i = 0; i < sd.nstloc(); i++ )
      {
        // divergence corrections
        double div_corr = 0.0;

        // SumExpG2 contribution
        const double  div_corr_1 = exfac * SumExpG2 * occ_ki_[i];
        div_corr += div_corr_1;
        const double e_div_corr_1 = -0.5 * div_corr_1 * occ_ki_[i];
        exchange_sum += e_div_corr_1;
        const double fac1 = 0.5 * exfac * occ_ki_[i] * occ_ki_[i];
        sigma_exhf_[0] += ( e_div_corr_1 + fac1 * sigma_sumexp[0] ) / omega;
        sigma_exhf_[1] += ( e_div_corr_1 + fac1 * sigma_sumexp[1] ) / omega;
        sigma_exhf_[2] += ( e_div_corr_1 + fac1 * sigma_sumexp[2] ) / omega;
        sigma_exhf_[3] += ( fac1 * sigma_sumexp[3] ) / omega;
        sigma_exhf_[4] += ( fac1 * sigma_sumexp[4] ) / omega;
        sigma_exhf_[5] += ( fac1 * sigma_sumexp[5] ) / omega;

        // rcut*rcut divergence correction
        if ( (vbasis_->mype() == 0) )
        {
          const double div_corr_2 = - alpha_sx_ * exfac *
             rcut_ * rcut_ * occ_ki_[i];
          div_corr += div_corr_2;
          const double e_div_corr_2 = -0.5 * div_corr_2 * occ_ki_[i];
          exchange_sum += e_div_corr_2;
          sigma_exhf_[0] += e_div_corr_2 / omega;
          sigma_exhf_[1] += e_div_corr_2 / omega;
          sigma_exhf_[2] += e_div_corr_2 / omega;
        }

        // analytical part
        if ( vbasis_->mype() == 0 )
        {
          const double div_corr_3 = - exfac * integ/vbz * occ_ki_[i];
          div_corr += div_corr_3;
          const double e_div_corr_3 = -0.5 * div_corr_3 * occ_ki_[i];
          exchange_sum += e_div_corr_3;
          // no contribution to stress
        }

        // contribution of divergence corrections to forces on wave functions
        if (dwf)
        {
          // sum the partial contributions to the correction for state i
          double tsum = div_corr;
          MPI_Allreduce(&tsum,&div_corr,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());

          // add correction to the derivatives of state i
          complex<double> *ps=c.valptr(i*c.mloc());
          complex<double> *pf=dc.valptr(i*dc.mloc());
          for ( int j = 0; j < dc.mloc(); j++ )
            pf[j] -= ps[j] * div_corr;
        }
      } // for i
    }

    // divergence corrections done

    if ( use_bisection_ )
    {
      bisection_[isp_loc]->backward(*uc_[isp_loc], *dwf->sd(isp_loc,0));
    }

  } // for isp_loc

  if ( use_bisection_ )
  {
    ostringstream ostr;
    for ( int ispin = 0; ispin < wfc_.nspin(); ++ispin )
    {
      int isp_loc = wfc_.isp_local(ispin);
      int isrc = -1;
      if ( isp_loc >= 0  )
      {
        if ( MPIdata::sd_rank() == 0 )
        {
          ostr.str("");
          isrc = MPIdata::rank();
          ostr << setprecision(10);
          ostr << " ExchangeOperator: bisection size: ispin=" << ispin
               << ": " << bisection_[isp_loc]->total_size() << endl;
          ostr << " ExchangeOperator: pair fraction:  ispin=" << ispin
               << ": " << bisection_[isp_loc]->pair_fraction() << endl;
        }
      }
      cout0(ostr.str(),isrc);
      MPI_Barrier(MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      cout << setprecision(3);
      cout << " ExchangeOperator: bisection time: "
           << tmb.real() << " s" << endl;
#if TIMING
      cout << setprecision(3);
      cout << " ExchangeOperator: bisection compute transform time: "
           << tmbtransf.real() << " s" << endl;
      cout << setprecision(3);
      cout << " ExchangeOperator: bisection compute localization time: "
           << tmbcomploc.real() << " s" << endl;
      cout << setprecision(3);
      cout << " ExchangeOperator: bisection forward time: "
           << tmbfwd.real() << " s" << endl;
#endif
    }
  }

  // sum contributions to the exchange energy
  double tsum = exchange_sum;
  MPI_Allreduce(&tsum,&exchange_sum,1,MPI_DOUBLE,MPI_SUM,MPIdata::comm());

  // accumulate stress tensor contributions
  if ( compute_stress )
  {
    double tsum[6];
    for ( int i = 0; i < 6; ++i )
      tsum[i] = sigma_exhf_[i];
    MPI_Allreduce(&tsum,&sigma_exhf_,1,MPI_DOUBLE,MPI_SUM,MPIdata::comm());
  }

  tm.stop();
  return exchange_sum;
}

////////////////////////////////////////////////////////////////////////////////
// Communication functions
////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::InitPermutation(void)
{
  // determine to whom we send and from whom we receive
  const int mycol = MPIdata::istb();
  const int npcol = MPIdata::nstb();

  colSendTo_ = ( mycol + 1 ) % npcol;
  colRecvFr_ = ( mycol - 1 + npcol ) % npcol;

  // Init communication for the number of states
  MPI_Send_init((void *) &nStatesKpi_, 1, MPI_INT, colSendTo_,
     NumberOfStates, MPIdata::st_comm(), &send_request_NumberOfStates_ );
  MPI_Recv_init((void *) &nNextStatesKpi_, 1, MPI_INT, colRecvFr_,
     NumberOfStates, MPIdata::st_comm(), &recv_request_NumberOfStates_ );
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::FreePermutation(void)
{
  // free permanent communications
  MPI_Request_free(&send_request_NumberOfStates_);
  MPI_Request_free(&recv_request_NumberOfStates_);
}

////////////////////////////////////////////////////////////////////////////////
// Permutation of the state numbers
// this function sets nNextStatesKpi_
void ExchangeOperator::SetNextPermutationStateNumber(void)
{
  // send the number of states to be send
  MPI_Start(&send_request_NumberOfStates_);
  MPI_Start(&recv_request_NumberOfStates_);

  // wait for the number of states to receive to be transmitted
  {
    MPI_Status status_recv;
    MPI_Status status_send;
    MPI_Wait( &recv_request_NumberOfStates_, &status_recv);
    MPI_Wait( &send_request_NumberOfStates_, &status_send);
  }
}

////////////////////////////////////////////////////////////////////////////////
// send the states in send_buf_states to the next column
// recieve the states from the previous column in state_kpi_
void ExchangeOperator::StartStatesPermutation(int mloc)
{
  // send the states
  if ( nStatesKpi_>0 )
  {
    wait_send_states_=1;
    MPI_Isend((void *) &send_buf_states_[0], 2*nStatesKpi_*mloc, MPI_DOUBLE,
      colSendTo_, States, MPIdata::st_comm(), &send_request_States_ );
  }
  else
  {
    wait_send_states_=0;
  }

  // receive the states
  if ( nNextStatesKpi_>0 )
  {
    wait_recv_states_=1;
    MPI_Irecv((void *) &state_kpi_[0], 2*nNextStatesKpi_*mloc, MPI_DOUBLE,
      colRecvFr_, States, MPIdata::st_comm(), &recv_request_States_ );
  }
  else
  {
    wait_recv_states_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteReceivingStates(int iRotationStep)
{
  // do something only if iRotationStep>0
  if ( iRotationStep != 0 && wait_recv_states_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &recv_request_States_, &status);

    // init reception flag
    wait_recv_states_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteSendingStates(int iRotationStep)
{
  // do something only if iRotationStep>0
  if ( iRotationStep != 0 && wait_send_states_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &send_request_States_, &status);

    // init reception flag
    wait_send_states_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
// StartForcesPermutation
// send the forces in send_buf_forces to the next column
// receive the forces from the previous column in force_kpi_
void ExchangeOperator::StartForcesPermutation(int mloc)
{
  // send the forces
  if ( nStatesKpi_>0 )
  {
    wait_send_forces_=1;
    MPI_Isend((void *) &send_buf_forces_[0], 2*nStatesKpi_*mloc, MPI_DOUBLE,
      colSendTo_, Forces, MPIdata::st_comm(), &send_request_Forces_ );
  }
  else
  {
    wait_send_forces_=0;
  }

  // receive the forces
  if ( nNextStatesKpi_>0 )
  {
    wait_recv_forces_=1;
    MPI_Irecv((void *) &force_kpi_[0], 2*nNextStatesKpi_*mloc, MPI_DOUBLE,
      colRecvFr_, Forces, MPIdata::st_comm(), &recv_request_Forces_ );
  }
  else
  {
    wait_recv_forces_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteReceivingForces(int iRotationStep)
{
  // do something only if iRotationStep>0
  if ( iRotationStep != 0 && wait_recv_forces_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &recv_request_Forces_, &status);

    // init reception flag
    wait_recv_forces_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteSendingForces(int iRotationStep)
{
  // do something only if iRotationStep>0
  if ( iRotationStep != 0 && wait_send_forces_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &send_request_Forces_, &status);

    // init reception flag
    wait_send_forces_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
// StartOccupationsPermutation
// store then send the occupations in send_buf_occupation_ to the next column
// receive the occupations from the previous column in occ_ki_
void ExchangeOperator::StartOccupationsPermutation(void)
{
  // store energies in the send buffer
  for ( int i=0; i<nStatesKpi_; i++ )
    send_buf_occupation_[i]=occ_ki_[i];

  // send the occupations
  if ( nStatesKpi_>0 )
  {
    wait_send_occupations_=1;
    MPI_Isend((void *) &send_buf_occupation_[0], nStatesKpi_, MPI_DOUBLE,
      colSendTo_, Occupation, MPIdata::st_comm(), &send_request_Occupation_ );
  }
  else
  {
    wait_send_occupations_=0;
  }

  // receive the occupations
  if ( nNextStatesKpi_>0 )
  {
    wait_recv_occupations_=1;
    MPI_Irecv((void *) &occ_ki_[0], nNextStatesKpi_, MPI_DOUBLE,
      colRecvFr_, Occupation, MPIdata::st_comm(), &recv_request_Occupation_ );
  }
  else
  {
    wait_recv_occupations_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteReceivingOccupations(int iRotationStep)
{
  if ( iRotationStep != 0 && wait_recv_occupations_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &recv_request_Occupation_, &status);

    // init reception flag
    wait_recv_occupations_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
void ExchangeOperator::CompleteSendingOccupations(int iRotationStep)
{
  if ( iRotationStep != 0 && wait_send_occupations_ )
  {
    // wait for the reception
    MPI_Status status;
    MPI_Wait( &send_request_Occupation_, &status);

    // init reception flag
    wait_send_occupations_=0;
  }
}

////////////////////////////////////////////////////////////////////////////////
// interaction potential
// vint(r) = alpha_sx_ * erf(mu_sx_*r)/r + beta_sx * erfc(mu_sx_*r)/r
//         = beta_sx_ / r + (alpha_sx_ - beta_sx_) * erf(mu_sx_*r)/r
// Fourier transform:
// vint(g2) = ( beta_sx + (alpha_sx-beta_sx) * exp( -g2 / 4*mu^2 ) ) / g2
double ExchangeOperator::vint(double g2)
{
  // this function should not be called with mu_sx_ == 0.0 if alpha != beta
  if ( alpha_sx_ != beta_sx_ ) assert(mu_sx_ != 0.0);
  // check if zero potential
  if ( ( alpha_sx_ == 0.0 ) && ( beta_sx_ == 0.0 ) )
    return 0.0;

  // Treat the Coulomb potential case separately (alpha_sx_ == beta_sx_)
  if ( alpha_sx_ == beta_sx_ )
  {
    if ( g2 == 0 )
      return 0.0;
    else
      return alpha_sx_ / g2;
  }
  else
  {
    const double fac = 0.25 / ( mu_sx_ * mu_sx_ );
    const double x = g2 * fac;
    if ( g2 == 0 )
      // return only the finite limit as g2 -> 0
      return - ( alpha_sx_ - beta_sx_ ) * fac;
    else if ( g2 < 1.e-6 )
      // Use Taylor expansion of the regular part near origin
      return alpha_sx_ / g2 + fac * beta_sx_ * ( 1.0 - 0.5 * x );
    else
      return ( beta_sx_ + ( alpha_sx_ - beta_sx_ ) * exp(-x) ) / g2;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Derivative of the interaction potential vint(g2) w.r.t g2
// dvint(g2) = d (vint(g2)) / d g2
double ExchangeOperator::dvint(double g2)
{
  // this function should not be called with mu_sx_ == 0.0 if alpha != beta
  if ( alpha_sx_ != beta_sx_ ) assert(mu_sx_ != 0.0);

  // check if zero potential
  if ( ( alpha_sx_ == 0.0 ) && ( beta_sx_ == 0.0 ) )
    return 0.0;

  // Treat the Coulomb potential case separately (alpha_sx_ == beta_sx_)
  if ( alpha_sx_ == beta_sx_ )
  {
    // Coulomb potential with prefactor alpha_sx (= beta_sx)
    if ( g2 == 0 )
      return 0.0;
    else
      return - alpha_sx_ / ( g2 * g2 );
  }
  else
  {
    const double fac = 0.25 / ( mu_sx_ * mu_sx_ );
    const double x = g2 * fac;
    const double third = 1.0 / 3.0;
    if ( g2 == 0 )
      // return finite part of the limit only for g2 -> 0
      return 0.5 * ( alpha_sx_ - beta_sx_ ) * fac * fac;
    else if ( g2 < 1e-6 )
      // Use Taylor expansion of regular term near origin
      // return beta_sx_ * ( -0.5 + fac * g2 * third ) * fac * fac;
      return - alpha_sx_ / ( g2 * g2 ) +
             0.5 * ( alpha_sx_ - beta_sx_ ) * fac * fac -
             ( alpha_sx_ - beta_sx_ ) * g2 * third * fac * fac * fac;
    else
      // exact derivative
      return - ( beta_sx_ / g2 +
        ( alpha_sx_ - beta_sx_ ) * exp(-x) * ( fac + 1.0/g2 ) ) / g2;
  }
}
