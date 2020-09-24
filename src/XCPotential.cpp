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
// XCPotential.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "XCPotential.h"
#include "LDAFunctional.h"
#include "VWNFunctional.h"
#include "PBEFunctional.h"
#include "BLYPFunctional.h"
#include "HSEFunctional.h"
#include "RSHFunctional.h"
#include "B3LYPFunctional.h"
#include "BHandHLYPFunctional.h"
#include "SCANFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "blas.h" // daxpy, dcopy
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
XCPotential::XCPotential(const ChargeDensity& cd, const string functional_name,
  const Sample& s): cd_(cd), vft_(*cd_.vft()), vbasis_(*cd_.vbasis()), s_(s)
{
  // copy arrays to resize rhototal_r_ and rhototal_g_
  rhototal_r_ = cd_.rhor;
  rhototal_g_ = cd_.rhog;
  if ( functional_name == "LDA" )
  {
    xcf_ = new LDAFunctional(rhototal_r_);
  }
  else if ( functional_name == "VWN" )
  {
    xcf_ = new VWNFunctional(rhototal_r_);
  }
  else if ( functional_name == "PBE" )
  {
    xcf_ = new PBEFunctional(rhototal_r_);
  }
  else if ( functional_name == "BLYP" )
  {
    xcf_ = new BLYPFunctional(rhototal_r_);
  }
  else if ( functional_name == "SCAN" )
  {
    xcf_ = new SCANFunctional(cd_.rhor);
  }
  else if ( functional_name == "PBE0" )
  {
    const double x_coeff = 1.0 - s_.ctrl.alpha_PBE0;
    const double c_coeff = 1.0;
    xcf_ = new PBEFunctional(rhototal_r_,x_coeff,c_coeff);
  }
  else if ( functional_name == "HSE" )
  {
    xcf_ = new HSEFunctional(rhototal_r_);
  }
  else if ( functional_name == "RSH" )
  {
    xcf_ = new RSHFunctional(rhototal_r_,s_.ctrl.alpha_RSH,s_.ctrl.beta_RSH,
           s_.ctrl.mu_RSH);
  }
  else if ( functional_name == "B3LYP" )
  {
    xcf_ = new B3LYPFunctional(rhototal_r_);
  }
  else if ( functional_name == "BHandHLYP" )
  {
    xcf_ = new BHandHLYPFunctional(rhototal_r_);
  }
  else
  {
    throw XCPotentialException("unknown functional name");
  }
  nspin_ = cd_.rhor.size();
  ngloc_ = vbasis_.localsize();
  np012loc_ = vft_.np012loc();

  if ( isGGA() )
  {
    vxctmp.resize(nspin_);
    for ( int ispin = 0; ispin < nspin_; ispin++ )
      vxctmp[ispin].resize(np012loc_);
    tmpr.resize(np012loc_);
  }
}

////////////////////////////////////////////////////////////////////////////////
XCPotential::~XCPotential(void)
{
  delete xcf_;
}

////////////////////////////////////////////////////////////////////////////////
bool XCPotential::isGGA(void)
{
  return xcf_->isGGA();
}

////////////////////////////////////////////////////////////////////////////////
bool XCPotential::isMeta(void)
{
  return xcf_->isMeta();
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::update(vector<vector<double> >& vr)
{
  // compute exchange-correlation energy and vxc potential vr[ispin][ir]

  // Input: total electronic density in:
  //   vector<vector<double> >           cd_.rhor[ispin][ir] (real space)
  //   vector<vector<complex<double> > > cd_.rhog[ispin][ig] (Fourier coeffs)
  // The array cd_.rhog is only used if xcf->isGGA() == true
  // to compute the density gradients

  // if a non-linear core correction is included,
  // rhototal_r_ = cd_.rhor+cd_.rhocore_r. Otherwise
  // rhototal_r_ = cd_.rhor

  // Output: (through member function xcf())
  //
  // exc_, dxc_
  //
  // LDA Functional:
  //   exc_, dxc_
  //   spin unpolarized: xcf()->exc, xcf()->vxc1
  //   spin polarized:   xcf()->exc, xcf()->vxc1_up, xcf()->vxc1_dn
  //
  // GGA Functional: (through member function xcf())
  //   exc_, dxc_
  //   spin unpolarized: xcf()->exc, xcf()->vxc1, xcf()->vxc2
  //   spin polarized:   xcf()->exc_up, xcf()->exc_dn,
  //                     xcf()->vxc1_up, xcf()->vxc1_dn
  //                     xcf()->vxc2_upup, xcf()->vxc2_dndn,
  //                     xcf()->vxc2_updn, xcf()->vxc2_dnup

  rhototal_r_ = cd_.rhor;
  rhototal_g_ = cd_.rhog;
  // test if a non-linear core correction is used
  // if so, add core density to rhototal_r_ and rhototal_g_
  if ( !cd_.rhocore_r.empty() )
  {
    // add core charge
    // note: if nspin==2, the cd_.rhocore_{rg} vectors each
    // contain half of the total core charge
    for ( int ispin = 0; ispin < rhototal_r_.size(); ispin++ )
    {
      //for ( int i = 0; i < rhototal_r_[ispin].size(); i++ )
      //  rhototal_r_[ispin][i] += cd_.rhocore_r[i];
      int len = rhototal_r_[ispin].size();
      int inc1 = 1;
      double one = 1.0;
      daxpy(&len,&one,(double*)&cd_.rhocore_r[0],&inc1,
            &rhototal_r_[ispin][0],&inc1);

      //for ( int i = 0; i < rhototal_g_[ispin].size(); i++ )
      //  rhototal_g_[ispin][i] += cd_.rhocore_g[i];
      len = 2*rhototal_g_[ispin].size();
      daxpy(&len,&one,(double*)&cd_.rhocore_g[0],&inc1,
            (double*)&rhototal_g_[ispin][0],&inc1);
    }
  }

  if ( !isGGA() )
  {
    // LDA functional

    xcf_->setxc();

    exc_ = 0.0;
    dxc_ = 0.0;
    const double *const e = xcf_->exc;
    const int size = xcf_->np();

    if ( nspin_ == 1 )
    {
      // unpolarized
      const double *const rh = xcf_->rho;
      const double *const v = xcf_->vxc1;
      for ( int i = 0; i < size; i++ )
      {
        const double e_i = e[i];
        const double v_i = v[i];
        const double rh_i = rh[i];
        exc_ += rh_i * e_i;
        dxc_ += rh_i * ( e_i - v_i );
        vr[0][i] = v_i;
      }
    }
    else
    {
      // spin polarized
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      const double *const v_up = xcf_->vxc1_up;
      const double *const v_dn = xcf_->vxc1_dn;
      for ( int i = 0; i < size; i++ )
      {
        const double r_i = rh_up[i] + rh_dn[i];
        exc_ += r_i * e[i];
        dxc_ += r_i * e[i] - rh_up[i] * v_up[i] - rh_dn[i] * v_dn[i];
        vr[0][i] = v_up[i];
        vr[1][i] = v_dn[i];
      }
    }
    double sum[2],tsum[2];
    sum[0] = exc_ * vbasis_.cell().volume() / vft_.np012();
    sum[1] = dxc_ * vbasis_.cell().volume() / vft_.np012();
    MPI_Allreduce(&sum,&tsum,2,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
    exc_ = tsum[0];
    dxc_ = tsum[1];
  }
  else
  {
    // GGA functional
    exc_ = 0.0;

    // compute grad_rho
    const double omega_inv = 1.0 / vbasis_.cell().volume();
    vector<complex<double> > tmp0(ngloc_), tmp1(ngloc_);
    if ( nspin_ == 1 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          /* i*G_j*c(G) */
          tmp0[ig] = complex<double>(0.0,omega_inv*gxj[ig])*rhototal_g_[0][ig];
        }
        vft_.backward(&tmp0[0],&tmpr[0]);
        int inc2=2, inc1=1;
        double *grj = xcf_->grad_rho[j];
        dcopy(&np012loc_,(double*)&tmpr[0],&inc2,grj,&inc1);
      }
    }
    else
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        const complex<double>* rhg0 = &rhototal_g_[0][0];
        const complex<double>* rhg1 = &rhototal_g_[1][0];
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          /* i*G_j*c(G) */
          const complex<double> igxj(0.0,omega_inv*gxj[ig]);
          const complex<double> c0 = *rhg0++;
          const complex<double> c1 = *rhg1++;
          tmp0[ig] = igxj * c0;
          tmp1[ig] = igxj * c1;
        }
        vft_.backward(&tmp0[0],&tmp1[0],&tmpr[0]);
        double *grj_up = xcf_->grad_rho_up[j];
        double *grj_dn = xcf_->grad_rho_dn[j];
        int inc2=2, inc1=1;
        double* p = (double*) &tmpr[0];
        dcopy(&np012loc_,p,  &inc2,grj_up,&inc1);
        dcopy(&np012loc_,p+1,&inc2,grj_dn,&inc1);
      } // j
    }
    if ( isMeta() )
    {
      // compute tau
      if (nspin_ == 1)
      {
        cd_.update_taur(xcf_->tau);
      }
      else
      {
        cd_.update_taur(xcf_->tau_up, xcf_->tau_dn);
      }
    }

    xcf_->setxc();

    // compute xc potential
    // take divergence of grad(rho)*vxc2

    // compute components of grad(rho) * vxc2
    if ( nspin_ == 1 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        const double *const grj = xcf_->grad_rho[j];
        const double *const v2 = xcf_->vxc2;
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          tmpr[ir] = grj[ir] * v2[ir];
        }
        // derivative
        vft_.forward(&tmpr[0],&tmp0[0]);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          tmp0[ig] *= complex<double>(0.0,gxj[ig]);
        }
        // back to real space
        vft_.backward(&tmp0[0],&tmpr[0]);
        // accumulate div(vxc2*grad_rho) in vxctmp
        double one = 1.0;
        int inc1 = 1, inc2 = 2;
        if ( j == 0 )
        {
          dcopy(&np012loc_,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
        }
        else
        {
          daxpy(&np012loc_,&one,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
        }
      }
    }
    else
    {
      double *v2_upup = xcf_->vxc2_upup;
      double *v2_updn = xcf_->vxc2_updn;
      double *v2_dnup = xcf_->vxc2_dnup;
      double *v2_dndn = xcf_->vxc2_dndn;
      for ( int j = 0; j < 3; j++ )
      {
        const double *gxj = vbasis_.gx_ptr(j);
        const double *grj_up = xcf_->grad_rho_up[j];
        const double *grj_dn = xcf_->grad_rho_dn[j];
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          const double re = v2_upup[ir] * grj_up[ir] + v2_updn[ir] * grj_dn[ir];
          const double im = v2_dnup[ir] * grj_up[ir] + v2_dndn[ir] * grj_dn[ir];
          tmpr[ir] = complex<double>(re,im);
        }
        // derivative
        vft_.forward(&tmpr[0],&tmp0[0],&tmp1[0]);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          const complex<double> igxj(0.0,gxj[ig]);
          tmp0[ig] *= igxj;
          tmp1[ig] *= igxj;
        }
        vft_.backward(&tmp0[0],&tmp1[0],&tmpr[0]);
        // accumulate div(vxc2*grad_rho) in vxctmp
        double one = 1.0;
        int inc1 = 1, inc2 = 2;
        double* p = (double*) &tmpr[0];
        if ( j == 0 )
        {
          dcopy(&np012loc_,p  ,&inc2,&vxctmp[0][0],&inc1);
          dcopy(&np012loc_,p+1,&inc2,&vxctmp[1][0],&inc1);
        }
        else
        {
          daxpy(&np012loc_,&one,p  ,&inc2,&vxctmp[0][0],&inc1);
          daxpy(&np012loc_,&one,p+1,&inc2,&vxctmp[1][0],&inc1);
        }
      } // j
    }

    // xc potential vr[i]
    // div(vxc2*grad_rho) is stored in vxctmp[ispin][ir]

    double esum=0.0;
    double dsum=0.0;
    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const rh = xcf_->rho;
      {
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          const double e_i = e[ir];
          const double rh_i = rh[ir];
          const double v_i = v1[ir] + vxctmp[0][ir];
          esum += rh_i * e_i;
          dsum += rh_i * ( e_i - v_i );
          vr[0][ir] = v_i;
        }
      }
    }
    else
    {
      const double *const v1_up = xcf_->vxc1_up;
      const double *const v1_dn = xcf_->vxc1_dn;
      const double *const eup = xcf_->exc_up;
      const double *const edn = xcf_->exc_dn;
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        const double r_up_i = rh_up[ir];
        const double r_dn_i = rh_dn[ir];
        esum += r_up_i * eup[ir] + r_dn_i * edn[ir];
        const double v_up = v1_up[ir] + vxctmp[0][ir];
        const double v_dn = v1_dn[ir] + vxctmp[1][ir];
        dsum += r_up_i * ( eup[ir] - v_up ) + r_dn_i * ( edn[ir] - v_dn );
        vr[0][ir] = v_up;
        vr[1][ir] = v_dn;
      }
    }
    double sum[2], tsum[2];
    sum[0] = esum * vbasis_.cell().volume() / vft_.np012();
    sum[1] = dsum * vbasis_.cell().volume() / vft_.np012();
    MPI_Allreduce(&sum,&tsum,2,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
    exc_ = tsum[0];
    dxc_ = tsum[1];
  }

  if ( isMeta() )
  {
    if ( nspin_ == 1 )
    {
      double sum = 0.0;
      const double *const v3 = xcf_->vxc3;
      const double *const tau = xcf_->tau;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        sum += tau[ir] * v3[ir];
      }
      sum *= vbasis_.cell().volume() / vft_.np012();
      double tsum = 0.0;
      MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
      dxc_ -= tsum;
    }
    else
    {
      double sum_up = 0.0;
      double sum_dn = 0.0;
      const double *const v3_up = xcf_->vxc3_up;
      const double *const v3_dn = xcf_->vxc3_dn;
      const double *const tau_up = xcf_->tau_up;
      const double *const tau_dn = xcf_->tau_dn;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        sum_up += tau_up[ir] * v3_up[ir];
        sum_dn += tau_dn[ir] * v3_dn[ir];
      }
      sum_up *= vbasis_.cell().volume() / vft_.np012();
      sum_dn *= vbasis_.cell().volume() / vft_.np012();
      double tsum_up = 0.0;
      double tsum_dn = 0.0;
      MPI_Allreduce(&sum_up,&tsum_up,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
      MPI_Allreduce(&sum_dn,&tsum_dn,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());
      dxc_ -= (tsum_up + tsum_dn);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::compute_stress(valarray<double>& sigma_exc)
{
  // compute exchange-correlation contributions to the stress tensor

  if ( !isGGA() && !isMeta())
  {
    // LDA functional

    double dsum = 0.0;
    const double *const e = xcf_->exc;
    const int size = xcf_->np();

    if ( nspin_ == 1 )
    {
      // unpolarized
      const double *const rh = xcf_->rho;
      const double *const v = xcf_->vxc1;
      for ( int i = 0; i < size; i++ )
      {
        dsum += rh[i] * (e[i] - v[i]);
      }
    }
    else
    {
      // spin polarized
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      const double *const v_up = xcf_->vxc1_up;
      const double *const v_dn = xcf_->vxc1_dn;
      for ( int i = 0; i < size; i++ )
      {
        const double rh = rh_up[i] + rh_dn[i];
        dsum += rh * e[i] - rh_up[i] * v_up[i] - rh_dn[i] * v_dn[i];
      }
    }
    const double fac = 1.0 / vft_.np012();
    double sum, tsum;
    // Next line: factor omega in volume element cancels 1/omega in
    // definition of sigma_exc
    sum = - fac * dsum;
    MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());

    // Note: contribution to sigma_exc is a multiple of the identity
    sigma_exc[0] = tsum;
    sigma_exc[1] = tsum;
    sigma_exc[2] = tsum;
    sigma_exc[3] = 0.0;
    sigma_exc[4] = 0.0;
    sigma_exc[5] = 0.0;
  }
  else if ( isGGA() && !isMeta())
  {
    // GGA functional

    double dsum=0.0,sum0=0.0,sum1=0.0,sum2=0.0,
           sum3=0.0,sum4=0.0,sum5=0.0;
    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const v2 = xcf_->vxc2;
      const double *const rh = xcf_->rho;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        dsum += rh[ir] * ( e[ir] - v1[ir] );
        const double grx = xcf_->grad_rho[0][ir];
        const double gry = xcf_->grad_rho[1][ir];
        const double grz = xcf_->grad_rho[2][ir];
        const double grx2 = grx * grx;
        const double gry2 = gry * gry;
        const double grz2 = grz * grz;
        const double grad2 = grx2 + gry2 + grz2;
        const double v2t = v2[ir];
        sum0 += ( grad2 + grx2 ) * v2t;
        sum1 += ( grad2 + gry2 ) * v2t;
        sum2 += ( grad2 + grz2 ) * v2t;
        sum3 += grx * gry * v2t;
        sum4 += gry * grz * v2t;
        sum5 += grx * grz * v2t;
      }
    }
    else
    {
      const double *const v1_up = xcf_->vxc1_up;
      const double *const v1_dn = xcf_->vxc1_dn;
      const double *const v2_upup = xcf_->vxc2_upup;
      const double *const v2_updn = xcf_->vxc2_updn;
      const double *const v2_dnup = xcf_->vxc2_dnup;
      const double *const v2_dndn = xcf_->vxc2_dndn;
      const double *const eup = xcf_->exc_up;
      const double *const edn = xcf_->exc_dn;
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        const double r_up = rh_up[ir];
        const double r_dn = rh_dn[ir];
        dsum += r_up * ( eup[ir] - v1_up[ir] ) +
                r_dn * ( edn[ir] - v1_dn[ir] );

        const double grx_up = xcf_->grad_rho_up[0][ir];
        const double gry_up = xcf_->grad_rho_up[1][ir];
        const double grz_up = xcf_->grad_rho_up[2][ir];
        const double grx2_up = grx_up * grx_up;
        const double gry2_up = gry_up * gry_up;
        const double grz2_up = grz_up * grz_up;
        const double grad2_up = grx2_up + gry2_up + grz2_up;

        const double grx_dn = xcf_->grad_rho_dn[0][ir];
        const double gry_dn = xcf_->grad_rho_dn[1][ir];
        const double grz_dn = xcf_->grad_rho_dn[2][ir];
        const double grx2_dn = grx_dn * grx_dn;
        const double gry2_dn = gry_dn * gry_dn;
        const double grz2_dn = grz_dn * grz_dn;
        const double grad2_dn = grx2_dn + gry2_dn + grz2_dn;

        const double grad_up_grad_dn = grx_up * grx_dn +
                                       gry_up * gry_dn +
                                       grz_up * grz_dn;

        const double v2_upup_ir = v2_upup[ir];
        const double v2_updn_ir = v2_updn[ir];
        const double v2_dnup_ir = v2_dnup[ir];
        const double v2_dndn_ir = v2_dndn[ir];

        sum0 += v2_upup_ir * ( grad2_up + grx2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grx_up * grx_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grx_dn * grx_up ) +
                v2_dndn_ir * ( grad2_dn + grx2_dn );

        sum1 += v2_upup_ir * ( grad2_up + gry2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + gry_up * gry_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + gry_dn * gry_up ) +
                v2_dndn_ir * ( grad2_dn + gry2_dn );

        sum2 += v2_upup_ir * ( grad2_up + grz2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grz_up * grz_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grz_dn * grz_up ) +
                v2_dndn_ir * ( grad2_dn + grz2_dn );

        sum3 += v2_upup_ir * grx_up * gry_up +
                v2_updn_ir * grx_up * gry_dn +
                v2_dnup_ir * grx_dn * gry_up +
                v2_dndn_ir * grx_dn * gry_dn;

        sum4 += v2_upup_ir * gry_up * grz_up +
                v2_updn_ir * gry_up * grz_dn +
                v2_dnup_ir * gry_dn * grz_up +
                v2_dndn_ir * gry_dn * grz_dn;

        sum5 += v2_upup_ir * grx_up * grz_up +
                v2_updn_ir * grx_up * grz_dn +
                v2_dnup_ir * grx_dn * grz_up +
                v2_dndn_ir * grx_dn * grz_dn;
      }
    }
    double fac = 1.0 / vft_.np012();
    double sum[6],tsum[6];
    // Next line: factor omega in volume element cancels 1/omega in
    // definition of sigma_exc
    sum[0] = - fac * ( dsum + sum0 );
    sum[1] = - fac * ( dsum + sum1 );
    sum[2] = - fac * ( dsum + sum2 );
    sum[3] = - fac * sum3;
    sum[4] = - fac * sum4;
    sum[5] = - fac * sum5;
    MPI_Allreduce(sum,tsum,6,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());

    sigma_exc[0] = tsum[0];
    sigma_exc[1] = tsum[1];
    sigma_exc[2] = tsum[2];
    sigma_exc[3] = tsum[3];
    sigma_exc[4] = tsum[4];
    sigma_exc[5] = tsum[5];
  }
  else
  {
    // MetaGGA functional
    double dsum=0.0,sum0=0.0,sum1=0.0,sum2=0.0,
           sum3=0.0,sum4=0.0,sum5=0.0;
    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const v2 = xcf_->vxc2;
      const double *const v3 = xcf_->vxc3;
      const double *const rh = xcf_->rho;
      const double *const tau = xcf_->tau;

      //V3 Stress -V3 * sum_n[grad_i rho_n * grad_j rho_n]
      const Wavefunction& wf0 = s_.wf;

      vector<complex<double> > tmprx(np012loc_),tmpry(np012loc_),
                               tmprz(np012loc_);
      vector<double> v3txx(np012loc_,0.0), v3txy(np012loc_,0.0),
                     v3txz(np012loc_,0.0), v3tyy(np012loc_,0.0),
                     v3tyz(np012loc_,0.0), v3tzz(np012loc_,0.0);
      for ( int isp_loc = 0; isp_loc < wf0.nsp_loc(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < wf0.nkp_loc(); ++ikp_loc )
        {
          double omega_inv = 1.0 / wf0.cell().volume();
          const SlaterDet* sd = wf0.sd(isp_loc,ikp_loc);
          if ( sd->basis().real() )
          {
            const int ngwloc = sd->basis().localsize();
            vector<complex<double> > tmp0(ngwloc), tmp1(ngwloc);
            const int mloc = sd->c().mloc();
            const complex<double>* p = sd->c().cvalptr();
            const double *const gxjx = sd->basis().gx_ptr(0);
            const double *const gxjy = sd->basis().gx_ptr(1);
            const double *const gxjz = sd->basis().gx_ptr(2);
            for ( int n = 0; n < sd->nstloc()-1; n++, n++ )
            {
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ0 = sd->occ(nn);
              double occ1 = sd->occ(nn+1);
              if ( ( occ0 + occ1 ) > 0.0 )
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmprx[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmpry[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight0 = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                       omega_inv * occ0;
                const double weight1 = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                       omega_inv * occ1;
                for ( int i = 0; i < np012loc_; i++ )
                {
                  v3txx[i] += weight0 * tmprx[i].real() * tmprx[i].real() +
                              weight1 * tmprx[i].imag() * tmprx[i].imag();
                  v3tyy[i] += weight0 * tmpry[i].real() * tmpry[i].real() +
                              weight1 * tmpry[i].imag() * tmpry[i].imag();
                  v3tzz[i] += weight0 * tmprz[i].real() * tmprz[i].real() +
                              weight1 * tmprz[i].imag() * tmprz[i].imag();
                  v3txy[i] += weight0 * tmprx[i].real() * tmpry[i].real() +
                              weight1 * tmprx[i].imag() * tmpry[i].imag();
                  v3txz[i] += weight0 * tmprx[i].real() * tmprz[i].real() +
                              weight1 * tmprx[i].imag() * tmprz[i].imag();
                  v3tyz[i] += weight0 * tmpry[i].real() * tmprz[i].real() +
                              weight1 * tmpry[i].imag() * tmprz[i].imag();
                }
              }
            }
            if ( sd->nstloc() % 2 != 0 )
            {
              const int n = sd->nstloc()-1;
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ = sd->occ(nn);
              if ( occ > 0.0 )
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprx[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpry[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                      omega_inv * occ;
                for ( int i = 0; i < np012loc_; i++ )
                {
                  v3txx[i] += weight * tmprx[i].real() * tmprx[i].real();
                  v3tyy[i] += weight * tmpry[i].real() * tmpry[i].real();
                  v3tzz[i] += weight * tmprz[i].real() * tmprz[i].real();
                  v3txy[i] += weight * tmprx[i].real() * tmpry[i].real();
                  v3txz[i] += weight * tmprx[i].real() * tmprz[i].real();
                  v3tyz[i] += weight * tmpry[i].real() * tmprz[i].real();
                }
              }
            }
          }
          else
          {
            const int ngwloc = sd->basis().localsize();
            vector<complex<double> > tmp0(ngwloc), tmp1(ngwloc);
            const int mloc = sd->c().mloc();
            const complex<double>* p = sd->c().cvalptr();
            const double *const kpgxjx = sd->basis().kpgx_ptr(0);
            const double *const kpgxjy = sd->basis().kpgx_ptr(1);
            const double *const kpgxjz = sd->basis().kpgx_ptr(2);
            for ( int n = 0; n < sd->nstloc(); n++)
            {
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ = sd->occ(nn);
              if (occ > 0.0)
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjx[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprx[0]);

                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjy[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpry[0]);

                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjz[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                      omega_inv * occ;
                for ( int i = 0; i < np012loc_; i++ )
                {
                  v3txx[i] += weight * real(conj(tmprx[i]) * tmprx[i]);
                  v3tyy[i] += weight * real(conj(tmpry[i]) * tmpry[i]);
                  v3tzz[i] += weight * real(conj(tmprz[i]) * tmprz[i]);
                  v3txy[i] += weight * real(conj(tmprx[i]) * tmpry[i]);
                  v3txz[i] += weight * real(conj(tmprx[i]) * tmprz[i]);
                  v3tyz[i] += weight * real(conj(tmpry[i]) * tmprz[i]);
                }
              }
            } // n
          } // real/complex
        } // ikp_loc
      } // isp_loc

      vector<double> v3tmp(np012loc_);
      v3tmp = v3txx;
      MPI_Allreduce(&v3tmp[0],&v3txx[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyy;
      MPI_Allreduce(&v3tmp[0],&v3tyy[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tzz;
      MPI_Allreduce(&v3tmp[0],&v3tzz[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txy;
      MPI_Allreduce(&v3tmp[0],&v3txy[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txz;
      MPI_Allreduce(&v3tmp[0],&v3txz[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyz;
      MPI_Allreduce(&v3tmp[0],&v3tyz[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());

      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        //V1 Stress rho * (Exc - V1)
        dsum += rh[ir] * (e[ir] - v1[ir] );
        //V2 Stress V2 * (delta_ij * grad^2 + grad_i * grad_j)
        const double grx = xcf_->grad_rho[0][ir];
        const double gry = xcf_->grad_rho[1][ir];
        const double grz = xcf_->grad_rho[2][ir];
        const double grx2 = grx * grx;
        const double gry2 = gry * gry;
        const double grz2 = grz * grz;
        const double grad2 = grx2 + gry2 + grz2;
        const double v2t = v2[ir];
        const double v3t = v3[ir];
        const double taur = tau[ir];
        sum0 += ( grad2 + grx2 ) * v2t - v3t * (taur + v3txx[ir]);
        sum1 += ( grad2 + gry2 ) * v2t - v3t * (taur + v3tyy[ir]);
        sum2 += ( grad2 + grz2 ) * v2t - v3t * (taur + v3tzz[ir]);
        sum3 += grx * gry * v2t - v3t * (v3txy[ir]);
        sum4 += gry * grz * v2t - v3t * (v3tyz[ir]);
        sum5 += grx * grz * v2t - v3t * (v3txz[ir]);
      }
    }
    else
    {
      const double *const v1_up = xcf_->vxc1_up;
      const double *const v1_dn = xcf_->vxc1_dn;
      const double *const v2_upup = xcf_->vxc2_upup;
      const double *const v2_updn = xcf_->vxc2_updn;
      const double *const v2_dnup = xcf_->vxc2_dnup;
      const double *const v2_dndn = xcf_->vxc2_dndn;
      const double *const v3_up = xcf_->vxc3_up;
      const double *const v3_dn = xcf_->vxc3_dn;
      const double *const eup = xcf_->exc_up;
      const double *const edn = xcf_->exc_dn;
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      const double *const tau_up = xcf_->tau_up;
      const double *const tau_dn = xcf_->tau_dn;

      //V3 Stress -V3 * sum_n[grad_i rho_n * grad_j rho_n]
      const Wavefunction& wf0 = s_.wf;

      vector<complex<double> > tmprx(np012loc_),tmpry(np012loc_),
                               tmprz(np012loc_);
      vector<double> v3txxup(np012loc_,0.0), v3txyup(np012loc_,0.0),
                     v3txzup(np012loc_,0.0), v3tyyup(np012loc_,0.0),
                     v3tyzup(np012loc_,0.0), v3tzzup(np012loc_,0.0),
                     v3txxdn(np012loc_,0.0), v3txydn(np012loc_,0.0),
                     v3txzdn(np012loc_,0.0), v3tyydn(np012loc_,0.0),
                     v3tyzdn(np012loc_,0.0), v3tzzdn(np012loc_,0.0);
      for ( int isp_loc = 0; isp_loc < wf0.nsp_loc(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < wf0.nkp_loc(); ++ikp_loc )
        {
          double omega_inv = 1.0 / wf0.cell().volume();
          const SlaterDet* sd = wf0.sd(isp_loc,ikp_loc);
          if ( sd->basis().real() )
          {
            const int ngwloc = sd->basis().localsize();
            vector<complex<double> > tmp0(ngwloc), tmp1(ngwloc);
            const int mloc = sd->c().mloc();
            const complex<double>* p = sd->c().cvalptr();
            const double *const gxjx = sd->basis().gx_ptr(0);
            const double *const gxjy = sd->basis().gx_ptr(1);
            const double *const gxjz = sd->basis().gx_ptr(2);
            for ( int n = 0; n < sd->nstloc()-1; n++, n++ )
            {
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ0 = sd->occ(nn);
              double occ1 = sd->occ(nn+1);
              if ( ( occ0 + occ1 ) > 0.0 )
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmprx[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmpry[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+n*mloc];
                  tmp1[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+(n+1)*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight0 = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                       omega_inv * occ0;
                const double weight1 = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                       omega_inv * occ1;
                if ( wf0.isp_global(isp_loc)==0 )
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxup[i] += weight0 * tmprx[i].real() * tmprx[i].real() +
                                  weight1 * tmprx[i].imag() * tmprx[i].imag();
                    v3tyyup[i] += weight0 * tmpry[i].real() * tmpry[i].real() +
                                  weight1 * tmpry[i].imag() * tmpry[i].imag();
                    v3tzzup[i] += weight0 * tmprz[i].real() * tmprz[i].real() +
                                  weight1 * tmprz[i].imag() * tmprz[i].imag();
                    v3txyup[i] += weight0 * tmprx[i].real() * tmpry[i].real() +
                                  weight1 * tmprx[i].imag() * tmpry[i].imag();
                    v3txzup[i] += weight0 * tmprx[i].real() * tmprz[i].real() +
                                  weight1 * tmprx[i].imag() * tmprz[i].imag();
                    v3tyzup[i] += weight0 * tmpry[i].real() * tmprz[i].real() +
                                  weight1 * tmpry[i].imag() * tmprz[i].imag();
                  }
                }
                else
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxdn[i] += weight0 * tmprx[i].real() * tmprx[i].real() +
                                  weight1 * tmprx[i].imag() * tmprx[i].imag();
                    v3tyydn[i] += weight0 * tmpry[i].real() * tmpry[i].real() +
                                  weight1 * tmpry[i].imag() * tmpry[i].imag();
                    v3tzzdn[i] += weight0 * tmprz[i].real() * tmprz[i].real() +
                                  weight1 * tmprz[i].imag() * tmprz[i].imag();
                    v3txydn[i] += weight0 * tmprx[i].real() * tmpry[i].real() +
                                  weight1 * tmprx[i].imag() * tmpry[i].imag();
                    v3txzdn[i] += weight0 * tmprx[i].real() * tmprz[i].real() +
                                  weight1 * tmprx[i].imag() * tmprz[i].imag();
                    v3tyzdn[i] += weight0 * tmpry[i].real() * tmprz[i].real() +
                                  weight1 * tmpry[i].imag() * tmprz[i].imag();
                  }
                }
              }
            }
            if ( sd->nstloc() % 2 != 0 )
            {
              const int n = sd->nstloc()-1;
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ = sd->occ(nn);
              if ( occ > 0.0 )
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjx[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprx[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjy[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpry[0]);
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,gxjz[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                      omega_inv * occ;
                if ( wf0.isp_global(isp_loc)==0 )
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxup[i] += weight * (tmprx[i].real() * tmprx[i].real());
                    v3tyyup[i] += weight * (tmpry[i].real() * tmpry[i].real());
                    v3tzzup[i] += weight * (tmprz[i].real() * tmprz[i].real());
                    v3txyup[i] += weight * (tmprx[i].real() * tmpry[i].real());
                    v3txzup[i] += weight * (tmprx[i].real() * tmprz[i].real());
                    v3tyzup[i] += weight * (tmpry[i].real() * tmprz[i].real());
                  }
                }
                else
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxdn[i] += weight * (tmprx[i].real() * tmprx[i].real());
                    v3tyydn[i] += weight * (tmpry[i].real() * tmpry[i].real());
                    v3tzzdn[i] += weight * (tmprz[i].real() * tmprz[i].real());
                    v3txydn[i] += weight * (tmprx[i].real() * tmpry[i].real());
                    v3txzdn[i] += weight * (tmprx[i].real() * tmprz[i].real());
                    v3tyzdn[i] += weight * (tmpry[i].real() * tmprz[i].real());
                  }
                }
              }
            }
          }
          else
          {
            const int ngwloc = sd->basis().localsize();
            vector<complex<double> > tmp0(ngwloc);
            const int mloc = sd->c().mloc();
            const complex<double>* p = sd->c().cvalptr();
            const double *const kpgxjx = sd->basis().kpgx_ptr(0);
            const double *const kpgxjy = sd->basis().kpgx_ptr(1);
            const double *const kpgxjz = sd->basis().kpgx_ptr(2);
            for ( int n = 0; n < sd->nstloc(); n++)
            {
              double nn = sd->context().mycol() * sd->c().nb() + n;
              double occ = sd->occ(nn);
              if (occ > 0.0)
              {
                // Compute Grad_j psi_n(ikp)
                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjx[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_x psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprx[0]);

                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjy[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_y psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpry[0]);

                for ( int ig = 0; ig < ngwloc; ig++ )
                {
                  /* i*G_j1*c(G) */
                  tmp0[ig] = complex<double>(0.0,kpgxjz[ig]) * p[ig+n*mloc];
                }
                // Compute Grad_z psi_n(ikp)
                cd_.ft(ikp_loc)->backward(&tmp0[0],&tmprz[0]);

                // Compute V3 * Grad_j1 psi_n(ikp) * Grad_j2 psi_n(ikp)
                const double weight = wf0.weight(wf0.ikp_global(ikp_loc)) *
                                      omega_inv * occ;
                if ( wf0.isp_global(isp_loc)==0 )
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxup[i] += weight * real(conj(tmprx[i]) * tmprx[i]);
                    v3tyyup[i] += weight * real(conj(tmpry[i]) * tmpry[i]);
                    v3tzzup[i] += weight * real(conj(tmprz[i]) * tmprz[i]);
                    v3txyup[i] += weight * real(conj(tmprx[i]) * tmpry[i]);
                    v3txzup[i] += weight * real(conj(tmprx[i]) * tmprz[i]);
                    v3tyzup[i] += weight * real(conj(tmpry[i]) * tmprz[i]);
                  }
                }
                else
                {
                  for ( int i = 0; i < np012loc_; i++ )
                  {
                    v3txxdn[i] += weight * real(conj(tmprx[i]) * tmprx[i]);
                    v3tyydn[i] += weight * real(conj(tmpry[i]) * tmpry[i]);
                    v3tzzdn[i] += weight * real(conj(tmprz[i]) * tmprz[i]);
                    v3txydn[i] += weight * real(conj(tmprx[i]) * tmpry[i]);
                    v3txzdn[i] += weight * real(conj(tmprx[i]) * tmprz[i]);
                    v3tyzdn[i] += weight * real(conj(tmpry[i]) * tmprz[i]);
                  }
                }
              }
            }
          }
        }
      }

      vector<double> v3tmp(np012loc_);
      v3tmp = v3txxup;
      MPI_Allreduce(&v3tmp[0],&v3txxup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyyup;
      MPI_Allreduce(&v3tmp[0],&v3tyyup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tzzup;
      MPI_Allreduce(&v3tmp[0],&v3tzzup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txyup;
      MPI_Allreduce(&v3tmp[0],&v3txyup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txzup;
      MPI_Allreduce(&v3tmp[0],&v3txzup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyzup;
      MPI_Allreduce(&v3tmp[0],&v3tyzup[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());

      v3tmp = v3txxdn;
      MPI_Allreduce(&v3tmp[0],&v3txxdn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyydn;
      MPI_Allreduce(&v3tmp[0],&v3tyydn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tzzdn;
      MPI_Allreduce(&v3tmp[0],&v3tzzdn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txydn;
      MPI_Allreduce(&v3tmp[0],&v3txydn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3txzdn;
      MPI_Allreduce(&v3tmp[0],&v3txzdn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());
      v3tmp = v3tyzdn;
      MPI_Allreduce(&v3tmp[0],&v3tyzdn[0],np012loc_,MPI_DOUBLE,
                    MPI_SUM,MPIdata::st_kp_sp_comm());

      for ( int ir = 0; ir < np012loc_; ir++ )
      {
        const double r_up = rh_up[ir];
        const double r_dn = rh_dn[ir];
        dsum += r_up * ( eup[ir] - v1_up[ir] ) +
                r_dn * ( edn[ir] - v1_dn[ir] );

        const double grx_up = xcf_->grad_rho_up[0][ir];
        const double gry_up = xcf_->grad_rho_up[1][ir];
        const double grz_up = xcf_->grad_rho_up[2][ir];
        const double grx2_up = grx_up * grx_up;
        const double gry2_up = gry_up * gry_up;
        const double grz2_up = grz_up * grz_up;
        const double grad2_up = grx2_up + gry2_up + grz2_up;

        const double grx_dn = xcf_->grad_rho_dn[0][ir];
        const double gry_dn = xcf_->grad_rho_dn[1][ir];
        const double grz_dn = xcf_->grad_rho_dn[2][ir];
        const double grx2_dn = grx_dn * grx_dn;
        const double gry2_dn = gry_dn * gry_dn;
        const double grz2_dn = grz_dn * grz_dn;
        const double grad2_dn = grx2_dn + gry2_dn + grz2_dn;

        const double grad_up_grad_dn = grx_up * grx_dn +
                                       gry_up * gry_dn +
                                       grz_up * grz_dn;

        const double v2_upup_ir = v2_upup[ir];
        const double v2_updn_ir = v2_updn[ir];
        const double v2_dnup_ir = v2_dnup[ir];
        const double v2_dndn_ir = v2_dndn[ir];
        const double v3tup = v3_up[ir];
        const double v3tdn = v3_dn[ir];
        const double taurup = tau_up[ir];
        const double taurdn = tau_dn[ir];

        sum0 += v2_upup_ir * ( grad2_up + grx2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grx_up * grx_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grx_dn * grx_up ) +
                v2_dndn_ir * ( grad2_dn + grx2_dn ) -
                v3tup * (taurup + v3txxup[ir]) - v3tdn * (taurdn + v3txxdn[ir]);

        sum1 += v2_upup_ir * ( grad2_up + gry2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + gry_up * gry_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + gry_dn * gry_up ) +
                v2_dndn_ir * ( grad2_dn + gry2_dn ) -
                v3tup * (taurup + v3tyyup[ir]) - v3tdn * (taurdn + v3tyydn[ir]);

        sum2 += v2_upup_ir * ( grad2_up + grz2_up ) +
                v2_updn_ir * ( grad_up_grad_dn + grz_up * grz_dn ) +
                v2_dnup_ir * ( grad_up_grad_dn + grz_dn * grz_up ) +
                v2_dndn_ir * ( grad2_dn + grz2_dn ) -
                v3tup * (taurup + v3tzzup[ir]) - v3tdn * (taurdn + v3tzzdn[ir]);

        sum3 += v2_upup_ir * grx_up * gry_up +
                v2_updn_ir * grx_up * gry_dn +
                v2_dnup_ir * grx_dn * gry_up +
                v2_dndn_ir * grx_dn * gry_dn -
                v3tup * (v3txyup[ir]) - v3tdn * (v3txydn[ir]);

        sum4 += v2_upup_ir * gry_up * grz_up +
                v2_updn_ir * gry_up * grz_dn +
                v2_dnup_ir * gry_dn * grz_up +
                v2_dndn_ir * gry_dn * grz_dn -
                v3tup * (v3tyzup[ir]) - v3tdn * (v3tyzdn[ir]);

        sum5 += v2_upup_ir * grx_up * grz_up +
                v2_updn_ir * grx_up * grz_dn +
                v2_dnup_ir * grx_dn * grz_up +
                v2_dndn_ir * grx_dn * grz_dn -
                v3tup * (v3txzup[ir]) - v3tdn * (v3txzdn[ir]);
      }
    }
    double fac = 1.0 / vft_.np012();
    double sum[6],tsum[6];
    // Next line: factor omega in volume element cancels 1/omega in
    // definition of sigma_exc
    sum[0] = - fac * ( dsum + sum0 );
    sum[1] = - fac * ( dsum + sum1 );
    sum[2] = - fac * ( dsum + sum2 );
    sum[3] = - fac * sum3;
    sum[4] = - fac * sum4;
    sum[5] = - fac * sum5;
    MPI_Allreduce(sum,tsum,6,MPI_DOUBLE,MPI_SUM,MPIdata::g_comm());

    sigma_exc[0] = tsum[0];
    sigma_exc[1] = tsum[1];
    sigma_exc[2] = tsum[2];
    sigma_exc[3] = tsum[3];
    sigma_exc[4] = tsum[4];
    sigma_exc[5] = tsum[5];
  }
}

////////////////////////////////////////////////////////////////////////////////
void XCPotential::apply_meta_operator(Wavefunction& dwf)
{
  const Wavefunction& wf0 = s_.wf;

  for ( int isp_loc = 0; isp_loc < wf0.nsp_loc(); ++isp_loc )
  {
    const double* pvxc3 = 0;
    const int ispg = wf0.isp_global(isp_loc);
    if ( wf0.nspin() == 1 )
    {
      pvxc3 = xcf_->vxc3;
    }
    else
    {
      if ( ispg == 0 )
        pvxc3 = xcf_->vxc3_up;
      else
        pvxc3 = xcf_->vxc3_dn;
    }

    for ( int ikp_loc = 0; ikp_loc < wf0.nkp_loc(); ++ikp_loc )
    {
      const SlaterDet* sd = wf0.sd(isp_loc,ikp_loc);
      if ( sd->basis().real() )
      {
        const int ngwloc = sd->basis().localsize();
        vector<complex<double> > tmp0(ngwloc), tmp1(ngwloc);
        const int mloc = sd->c().mloc();
        const complex<double>* p = sd->c().cvalptr();
        complex<double>* dp = dwf.sd(isp_loc,ikp_loc)->c().valptr();
        for ( int n = 0; n < sd->nstloc()-1; n++, n++ )
        {
          for ( int j = 0; j < 3; j++ )
          {
            // Compute Grad_j psi_n(ikp)
            const double *const gxj = sd->basis().gx_ptr(j);
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              /* i*G_j*c(G) */
              tmp0[ig] = complex<double>(0.0,gxj[ig]) * p[ig+n*mloc];
              tmp1[ig] = complex<double>(0.0,gxj[ig]) * p[ig+(n+1)*mloc];
            }
            cd_.ft(ikp_loc)->backward(&tmp0[0],&tmp1[0],&tmpr[0]);
            // Compute V3 * Grad_j psi_n(ikp)
            for ( int i = 0; i < np012loc_; i++ )
              tmpr[i] *= pvxc3[i];
            // Transform to k-space
            cd_.ft(ikp_loc)->forward(&tmpr[0],&tmp0[0],&tmp1[0]);
            // Compute Div_j[V3 * Grad_j psi_n(ikp)]
            // Note -1/2 comes from definition of V3
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              /* i*G_j*c(G) */
              dp[ig+n*mloc] += -0.5 * complex<double>(0.0,gxj[ig]) * tmp0[ig];
              dp[ig+(n+1)*mloc] += -0.5 * complex<double>(0.0,gxj[ig]) *
                                   tmp1[ig];
            }
          }
        }
        if ( sd->nstloc() % 2 != 0 )
        {
          const int n = sd->nstloc()-1;
          for ( int j = 0; j < 3; j++ )
          {
            const double *const gxj = sd->basis().gx_ptr(j);
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              tmp0[ig] = complex<double>(0.0,gxj[ig]) * p[ig+n*mloc];
            }
            cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpr[0]);
            for ( int i = 0; i < np012loc_; i++ )
              tmpr[i] *= pvxc3[i];
            cd_.ft(ikp_loc)->forward(&tmpr[0],&tmp0[0]);
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              dp[ig+n*mloc] += -0.5 * complex<double>(0.0,gxj[ig]) * tmp0[ig];
            }
          }
        }
      }
      else
      {
        const int ngwloc = sd->basis().localsize();
        vector<complex<double> > tmp0(ngwloc), tmp1(ngwloc);
        const int mloc = sd->c().mloc();
        const complex<double>* p = sd->c().cvalptr();
        complex<double>* dp = dwf.sd(isp_loc,ikp_loc)->c().valptr();
        for ( int n = 0; n < sd->nstloc(); n++ )
        {
          for ( int j = 0; j < 3; j++ )
          {
            // Compute Grad_j psi_n(ikp)
            const double *const kpgxj = sd->basis().kpgx_ptr(j);
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              // i*(k+G)_j*c(G)
              tmp0[ig] = complex<double>(0.0,kpgxj[ig]) * p[ig+n*mloc];
            }
            cd_.ft(ikp_loc)->backward(&tmp0[0],&tmpr[0]);
            // Compute V3 * Grad_j psi_n(ikp)
            for ( int i = 0; i < np012loc_; i++ )
              tmpr[i] *= pvxc3[i];
            // Transform to k-space
            cd_.ft(ikp_loc)->forward(&tmpr[0],&tmp0[0]);
            // Compute Div_j[V3 * Grad_j psi_n(ikp)]
            // Note -1/2 comes from definition of V3
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              // i*(k+G)_j*c(G)
              dp[ig+n*mloc] += -0.5 * complex<double>(0.0,kpgxj[ig]) *
                               tmp0[ig];
            }
          }
        }
      }
    } // ikp_loc
  } // isp_loc
}
