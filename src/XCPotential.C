////////////////////////////////////////////////////////////////////////////////
//
// XCPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: XCPotential.C,v 1.2 2003-05-16 16:14:00 fgygi Exp $

#include "XCPotential.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "blas.h" // daxpy, dcopy
#include <cassert>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
XCPotential::XCPotential(ChargeDensity& cd, const string functional_name) :
cd_(cd), ctxt_(*cd.vcontext()), vft_(*cd_.vft()), vbasis_(*cd_.vbasis())
{
  if ( functional_name == "LDA" )
  {
    xcf_ = new LDAFunctional(cd_.rhor);
  }
  else if ( functional_name == "PBE" )
  {
    xcf_ = new PBEFunctional(cd_.rhor);
  }
  else if ( functional_name == "BLYP" )
  {
    xcf_ = new BLYPFunctional(cd_.rhor);
  }
  else
  {
    throw XCPotentialException("unknown functional name");
  }
  nspin_ = cd_.rhor.size();
  ngloc_ = vbasis_.localsize();
  np012loc_ = vft_.np012loc();
  
  if ( xcf_->isGGA() )
  {
    tmp1.resize(ngloc_);
    if ( nspin_ > 1 )
      tmp2.resize(ngloc_);
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
void XCPotential::update(vector<vector<double> >& vr, bool compute_stress)
{
  // compute exchange-correlation energy and add vxc potential to vr[ispin][ir]
  
  // Input: total electronic density in:
  //   vector<vector<double> >           cd_.rhor[ispin][ir] (real space)
  //   vector<vector<complex<double> > > cd_.rhog[ispin][ig] (Fourier coeffs)
  // The array cd_.rhog is only used if xcf->isGGA() == true
  // to compute the density gradients
  
  // Output: (through member function xcf())
  //
  // exc_, dxc, dxc0_, dxc1_, dxc2_
  //
  // LDA Functional:
  //   exc_, dxc
  //   spin unpolarized: xcf()->exc, xcf()->vxc1
  //   spin polarized:   xcf()->exc, xcf()->vxc1_up, xcf()->vxc1_dn
  //
  // GGA Functional: (through member function xcf())
  //   exc_, dxc, dxc0_, dxc1_, dxc2_
  //   spin unpolarized: xcf()->exc, xcf()->vxc1, xcf()->vxc2
  //   spin polarized:   xcf()->exc_up, xcf()->exc_dn, 
  //                     xcf()->vxc1_up, xcf()->vxc1_dn
  //                     xcf()->vxc2_upup, xcf()->vxc2_dndn, 
  //                     xcf()->vxc2_updn, xcf()->vxc2_dnup
  
  if ( !xcf_->isGGA() )
  {
    // LDA functional
 
    xcf_->setxc();
 
    exc_ = 0.0;
    dxc_ = 0.0;
    const double *const e = xcf_->exc;
 
    if ( nspin_ == 1 )
    {
      // unpolarized
      const double *const rh = xcf_->rho;
      const double *const v = xcf_->vxc1;
      const int size = xcf_->np();
      for ( int i = 0; i < size; i++ )
      {
        const double r = rh[i];
        const double tmp = r * e[i];
        exc_ += tmp;
        dxc_ += tmp - r * v[i];
        vr[0][i] += v[i];
      }
    }
    else
    {
      // spin polarized
      const double *const rh_up = xcf_->rho_up;
      const double *const rh_dn = xcf_->rho_dn;
      const double *const v_up = xcf_->vxc1_up;
      const double *const v_dn = xcf_->vxc1_dn;
      for ( int i = 0; i < np012loc_; i++ )                          
      {                                                         
        const double rh = rh_up[i] + rh_dn[i];                  
        const double tmp = rh * e[i];                           
        exc_ += tmp;                                            
        dxc_ += tmp - rh_up[i] * v_up[i] - rh_dn[i] * v_dn[i];  
        vr[0][i] += v_up[i];
        vr[1][i] += v_dn[i];
      }                                                         
    }
    double fac = vbasis_.cell().volume() / vft_.np012();
    exc_ *= fac;
    ctxt_.dsum(1,1,&exc_,1);
    //dexcda[0] = fac * dxc / al[0];
    //dexcda[1] = fac * dxc / al[1];
    //dexcda[2] = fac * dxc / al[2];
  }
  else
  {
    // GGA functional
    exc_ = 0.0;
    dxc_ = dxc0_ = dxc1_ = dxc2_ = 0.0;
    int size = xcf_->np();
    
    // compute grad_rho
    const double omega_inv = 1.0 / vbasis_.cell().volume();
    if ( nspin_ == 1 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          /* i*G_j*c(G) */
          tmp1[ig] = complex<double>(0.0,omega_inv*gxj[ig]) * cd_.rhog[0][ig];
        }
        vft_.backward(&tmp1[0],&tmpr[0]);
        int inc2=2, inc1=1;
        double *grj = xcf_->grad_rho[j];
        dcopy_(&np012loc_,(double*)&tmpr[0],&inc2,grj,&inc1);
      }
    }
    else
    {
      for ( int j = 0; j < 3; j++ )
      {
        const double *const gxj = vbasis_.gx_ptr(j);
        const complex<double>* rhg0 = &cd_.rhog[0][0];
        const complex<double>* rhg1 = &cd_.rhog[1][0];
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          /* i*G_j*c(G) */
          const complex<double> igxj(0.0,omega_inv*gxj[ig]);
          const complex<double> c0 = *rhg0++;
          const complex<double> c1 = *rhg1++;
          tmp1[ig] = igxj * c0; 
          tmp2[ig] = igxj * c1; 
        }
        vft_.backward(&tmp1[0],&tmp2[0],&tmpr[0]);
        double *grj_up = xcf_->grad_rho_up[j];
        double *grj_dn = xcf_->grad_rho_dn[j];
        int inc2=2, inc1=1;
        double* p = (double*) &tmpr[0];
        dcopy_(&np012loc_,p,  &inc2,grj_up,&inc1);
        dcopy_(&np012loc_,p+1,&inc2,grj_dn,&inc1);
      } // j
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
        vft_.forward(&tmpr[0],&tmp1[0]);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          tmp1[ig] *= complex<double>(0.0,gxj[ig]);
        }
        // back to real space
        vft_.backward(&tmp1[0],&tmpr[0]);
        // accumulate div(vxc2*grad_rho) in vxctmp
        double one = 1.0;
        int inc1 = 1, inc2 = 2;
        if ( j == 0 )
        {
          dcopy_(&np012loc_,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
        }
        else
        {
          daxpy_(&np012loc_,&one,(double*)&tmpr[0],&inc2,&vxctmp[0][0],&inc1);
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
        vft_.forward(&tmpr[0],&tmp1[0],&tmp2[0]);
        for ( int ig = 0; ig < ngloc_; ig++ )
        {
          // i*G_j*c(G)
          const complex<double> igxj(0.0,gxj[ig]);
          tmp1[ig] *= igxj;
          tmp2[ig] *= igxj;
        }
        vft_.backward(&tmp1[0],&tmp2[0],&tmpr[0]);
        // accumulate div(vxc2*grad_rho) in vxctmp
        double one = 1.0;
        int inc1 = 1, inc2 = 2;
        double* p = (double*) &tmpr[0];
        if ( j == 0 )
        {
          dcopy_(&np012loc_,p  ,&inc2,&vxctmp[0][0],&inc1);
          dcopy_(&np012loc_,p+1,&inc2,&vxctmp[1][0],&inc1);
        }
        else
        {
          daxpy_(&np012loc_,&one,p  ,&inc2,&vxctmp[0][0],&inc1);
          daxpy_(&np012loc_,&one,p+1,&inc2,&vxctmp[1][0],&inc1);
        }
      } // j
    }
    
    // add xc potential to local potential in vr[i]
    // div(vxc2*grad_rho) is stored in vxctmp[ispin][ir]

    if ( nspin_ == 1 )
    {
      const double *const e = xcf_->exc;
      const double *const v1 = xcf_->vxc1;
      const double *const v2 = xcf_->vxc2;
      const double *const rh = xcf_->rho;
      double esum=0.0,dsum=0.0,sum0=0.0,sum1=0.0,sum2=0.0;
      {
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          esum += rh[ir] * e[ir];
          vr[0][ir] += v1[ir] + vxctmp[0][ir];
          if ( compute_stress )
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
          }
        }
      }
      exc_ += esum;
      dxc_ += dsum;
      dxc0_ += sum0;
      dxc1_ += sum1;
      dxc2_ += sum2;
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
      double esum=0.0,dsum=0.0,sum0=0.0,sum1=0.0,sum2=0.0;
      {
        for ( int ir = 0; ir < np012loc_; ir++ )
        {
          const double r_up = rh_up[ir];
          const double r_dn = rh_dn[ir];
          esum += r_up * eup[ir] + r_dn * edn[ir];
          vr[0][ir] += v1_up[ir] + vxctmp[0][ir];
          vr[1][ir] += v1_dn[ir] + vxctmp[1][ir];
          if ( compute_stress )
          {
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

            sum0 += v2_upup[ir] * ( grad2_up + grx2_up ) +
            v2_updn[ir] * ( grad_up_grad_dn + grx_up * grx_dn ) +
            v2_dnup[ir] * ( grad_up_grad_dn + grx_up * grx_dn ) +
            v2_dndn[ir] * ( grad2_dn + grx2_dn );
 
            sum1 += v2_upup[ir] * ( grad2_up + gry2_up ) +
            v2_updn[ir] * ( grad_up_grad_dn + gry_up * gry_dn ) +
            v2_dnup[ir] * ( grad_up_grad_dn + gry_up * gry_dn ) +
            v2_dndn[ir] * ( grad2_dn + gry2_dn );
 
            sum2 += v2_upup[ir] * ( grad2_up + grz2_up ) +
            v2_updn[ir] * ( grad_up_grad_dn + grz_up * grz_dn ) +
            v2_dnup[ir] * ( grad_up_grad_dn + grz_up * grz_dn ) +
            v2_dndn[ir] * ( grad2_dn + grz2_dn );
 
          }
        }
      }
      exc_ += esum;
      dxc_ += dsum;
      dxc0_ += sum0;
      dxc1_ += sum1;
      dxc2_ += sum2;
    }
    double fac = vbasis_.cell().volume() / vft_.np012();
    exc_ *= fac;
    ctxt_.dsum(1,1,&exc_,1);
    //dexcda[0] = fac * (dxc + dxc0) / al[0];
    //dexcda[1] = fac * (dxc + dxc1) / al[1];
    //dexcda[2] = fac * (dxc + dxc2) / al[2];
    
  }
}
