////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.C,v 1.15 2003-10-02 17:40:34 fgygi Exp $

#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "StructureFactor.h"
#include "XCPotential.h"
#include "NonLocalPotential.h"

#include "Timer.h"
#include "blas.h"

#include <iostream>
#include <iomanip>
#include <algorithm> // fill()
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::EnergyFunctional(const Sample& s) : s_(s), cd_(s.wf)
{
  const AtomSet& atoms = s_.atoms;
  const Wavefunction& wf = s_.wf;
  const UnitCell& cell(wf.cell());
  double omega = cell.volume();
  double ecutv = 4 * wf.ecut();
  
  vbasis_ = cd_.vbasis();
  //cout << vbasis_->context().mype() << ": vbasis_->context() = " 
  //     << vbasis_->context() << endl;
  
  // define FT's on vbasis contexts
  
  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
  vft = cd_.vft();
  
  v_r.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    v_r[ispin].resize(vft->np012loc());
  }
  tmp_r.resize(vft->np012loc());

  if ( s_.ctxt_.onpe0() )
  {
    cout << "  <!-- EnergyFunctional: np0v,np1v,np2v: " << np0v << " "
         << np1v << " " << np2v << " -->" << endl;
    cout << "  <!-- EnergyFunctional: vft->np012(): "
         << vft->np012() << " -->" << endl;
  }
  
  int ngloc = vbasis_->localsize();
  //cout << " EnergyFunctional: ngloc: " << ngloc << endl;
  
  nsp_ = atoms.nsp();
  
  vps.resize(nsp_);
  dvps.resize(nsp_);
  rhops.resize(nsp_);
  
  zv_.resize(nsp_);
  rcps_.resize(nsp_);
  na_.resize(nsp_);
  namax_ = 0;
  
  for ( int is = 0; is < nsp_; is++ )
  {
    vps[is].resize(ngloc);
    dvps[is].resize(ngloc);
    rhops[is].resize(ngloc);
    if ( atoms.na(is) > namax_ ) namax_ = atoms.na(is);
  }
  
  xcp = new XCPotential(cd_,s_.ctrl.xc);
  nlp = new NonLocalPotential(s_.atoms, *wf.sd(0,0));
  
  // initialize Fourier coeffs of the potential 
  init();
  
  // FT for interpolation of wavefunctions on the fine grid
  ft.resize(wf.nkp());
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    ft[ikp] = cd_.ft(ikp);
  }
  
}

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::~EnergyFunctional(void)
{
  delete xcp;
  delete nlp;
  
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double EnergyFunctional::energy(bool compute_hpsi, Wavefunction& dwf,
              bool compute_forces, vector<vector<double> >& fion,
              bool compute_stress, UnitCell& dcell)
{
  const double fpi = 4.0 * M_PI;
  const double *g2i = vbasis_->g2i_ptr();

  const Wavefunction& wf = s_.wf;
  const UnitCell& cell(wf.cell());
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const int ngloc = vbasis_->localsize();
  vector<double> temp(ngloc);
  assert(wf.nspin()==1);
  
  if ( compute_forces )
  {
    for ( int is = 0; is < fion.size(); is++ )
    {
      for ( int i = 0; i < fion[is].size(); i++ )
      {
        fion[is][i] = 0.0;
      }
    }
  }
  
  if ( compute_hpsi )
  {
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        dwf.sd(ispin,ikp)->c().clear();
      }
    }
  }
        
  
  // kinetic energy
  tmap["ekin"].start();
  ekin_ = 0.0;
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      ekin_ += wf.sd(ispin,ikp)->ekin() * wf.weight(ikp);
    }
  }
  tmap["ekin"].stop();
  
  tmap["density"].start();
  cd_.update_density();
  tmap["density"].stop();
  
  // total electronic density: rhoelg = rho_up + rho_dn
  if ( wf.nspin() == 1 )
  {
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhoelg[ig] = omega_inv * cd_.rhog[0][ig];
    }
  }
  else
  {
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhoelg[ig] = omega_inv * ( cd_.rhog[0][ig] + cd_.rhog[1][ig] );
    }
  }
  
  // potential energy: integral of electronic charge times ionic local pot.
  int len=2*ngloc,inc1=1;
  eps_ = 2.0 * ddot_(&len,(double*)&rhoelg[0],&inc1,
         (double*)&vion_local_g[0],&inc1);
  // remove double counting for G=0
  if ( vbasis_->context().myrow() == 0 )
  {
    eps_ -= real(conj(rhoelg[0])*vion_local_g[0]);
  }
  eps_ *= omega;
  vbasis_->context().dsum(1,1,&eps_,1);

  // Hartree energy
  ehart_ = 0.0;
  for ( int ig = 0; ig < ngloc; ig++ )
  {
    const complex<double> tmp = rhoelg[ig] + rhopst[ig];
    ehart_ += norm(tmp) * g2i[ig];
    rhogt[ig] = tmp;
  }
  // factor 1/2 from definition of Ehart cancels with half sum over G
  ehart_ *= omega * 4.0 * M_PI;
  vbasis_->context().dsum(1,1,&ehart_,1);
  
  tmap["exc"].start();
  // XC energy and potential
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    fill(v_r[ispin].begin(),v_r[ispin].end(),0.0);
  // boolean in next line: compute_stress
  xcp->update(v_r,false);
  exc_ = xcp->exc();
  tmap["exc"].stop();
  
  // Non local energy
  // Note: next line for nspin==0, nkp==0 only
  tmap["nonlocal"].start();
  enl_ = nlp->energy(compute_hpsi,*dwf.sd(0,0),
    compute_forces, fion, compute_stress, dcell);
  tmap["nonlocal"].stop();

  ecoul_ = ehart_ + esr_ - eself_;
  etotal_ = ekin_ + eps_ + enl_ + ecoul_ + exc_;
  
  if ( compute_hpsi )
  {
  
    // compute vlocal_g = vion_local_g + vhart_g
    // where vhart_g = 4 * pi * (rhoelg + rhopst) * g2i
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      vlocal_g[ig] = vion_local_g[ig] + fpi * rhogt[ig] * g2i[ig];
    }
 
    // FT to tmpr_r
    vft->backward(&vlocal_g[0],&tmp_r[0]);
 
    // add local potential in tmp_r to v_r[ispin][i]
    // v_r contains the xc potential
    const int size = tmp_r.size();
    if ( wf.nspin() == 1 )
    {
      for ( int i = 0; i < size; i++ )
      {
        v_r[0][i] += real(tmp_r[i]);
      }
    }
    else
    {
      for ( int i = 0; i < size; i++ )
      {
        const double vloc = real(tmp_r[i]);
        v_r[0][i] += vloc;
        v_r[1][i] += vloc;
      }
    }
 
    tmap["hpsi"].start();
    assert(wf.nspin()==1);
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        if ( wf.sd(ispin,ikp) != 0 )
        {
          const SlaterDet& sd = *(wf.sd(ispin,ikp));
          SlaterDet& sdp = *(dwf.sd(ispin,ikp));
          const ComplexMatrix& c = sd.c();
          const Basis& wfbasis = sd.basis();
          ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
          const int mloc = cp.mloc();
          
          for ( int n = 0; n < sd.nstloc(); n++ )
          {
            // Laplacian
            for ( int ig = 0; ig < wfbasis.localsize(); ig++ )
            {
              cp[ig+mloc*n] += 0.5 * wfbasis.kpg2(ig) * c[ig+mloc*n];
            }
          }
          sd.rs_mul_add(*ft[ikp], &v_r[ispin][0], sdp);
        } // if sd(ispin,ikp) != 0
      }
    }
    tmap["hpsi"].stop();
  } // if compute_hpsi
  
  if ( compute_forces )
  {
    const int* idx = vbasis_->idx_ptr();
    const double* gx0 = vbasis_->gx_ptr(0);
    const double* gx1 = vbasis_->gx_ptr(1);
    const double* gx2 = vbasis_->gx_ptr(2);
    
    for ( int is = 0; is < nsp_; is++ )
    {
      for ( int ig = 0; ig < ngloc; ig++ )
      {
        double tmp = fpi * rhops[is][ig] * g2i[ig];
        vtemp[ig] =  tmp * conj(rhogt[ig]) + vps[is][ig] * conj(rhoelg[ig]);
      }
      memset((void*)&ftmp[0],0,3*namax_*sizeof(double));
      // loop over atoms of species is
      for ( int ia = 0; ia < na_[is]; ia++ )
      {
        double sum0=0.0,sum1=0.0,sum2=0.0;
        double *c0 = sf.cos0_ptr(is,ia);
        double *c1 = sf.cos1_ptr(is,ia);
        double *c2 = sf.cos2_ptr(is,ia);
        double *s0 = sf.sin0_ptr(is,ia);
        double *s1 = sf.sin1_ptr(is,ia);
        double *s2 = sf.sin2_ptr(is,ia);
        for ( int ig = 0; ig < ngloc; ig++ )
        {
          // compute Exp[igr] in 3D as a product of 1D contributions
//           complex<double> teigr = ei0[kv[3*ig+0]] *
//                                   ei1[kv[3*ig+1]] *
//                                   ei2[kv[3*ig+2]];
          const int iii = ig+ig+ig;
          const int kx = idx[iii];
          const int ky = idx[iii+1];
          const int kz = idx[iii+2];
          
          const double cos_a = c0[kx];
          const double cos_b = c1[ky];
          const double cos_c = c2[kz];
 
          const double sin_a = s0[kx];
          const double sin_b = s1[ky];
          const double sin_c = s2[kz];
 
          // Next line: exp(-i*gr) =
          // (cos_a - I sin_a)*(cos_b - I sin_b)*(cos_c - I sin_c)
          double teigr_re = 
            cos_a*cos_b*cos_c - sin_a*sin_b*cos_c -
            sin_a*cos_b*sin_c - cos_a*sin_b*sin_c;
          double teigr_im = 
            sin_a*sin_b*sin_c - sin_a*cos_b*cos_c -
            cos_a*sin_b*cos_c - cos_a*cos_b*sin_c;
                   
          /* fion is real */
          double tmp = teigr_re * vtemp[ig].imag() + 
                       teigr_im * vtemp[ig].real();

          sum0 += tmp * gx0[ig];
          sum1 += tmp * gx1[ig];
          sum2 += tmp * gx2[ig];
        }
        ftmp[3*ia]   = sum0;
        ftmp[3*ia+1] = sum1;
        ftmp[3*ia+2] = sum2;

      }
      
      int len = 3*na_[is];
      vbasis_->context().dsum(len,1,&ftmp[0],len);

      double fac = -2.0 * omega;
      for ( int i = 0; i < 3*na_[is]; i++ )
      {
        fion[is][i] += fion_esr[is][i] + fac * ftmp[i];
      }
    }
  }
  
  return etotal_;
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::init(void)
{
  int ngloc = vbasis_->localsize();
  vion_local_g.resize(ngloc);
  vlocal_g.resize(ngloc);
  vtemp.resize(ngloc);
  rhoelg.resize(ngloc);
  rhogt.resize(ngloc);
  rhopst.resize(ngloc);
  
  const Wavefunction& wf = s_.wf;
  const UnitCell& cell = wf.cell();
  double omega = cell.volume();
  assert(omega != 0.0);
  double omega_inv = 1.0 / omega;  
  
  const AtomSet& atoms = s_.atoms;
  
  const double *g;
  
  tau0.resize(nsp_);
  taum.resize(nsp_);
  fion_esr.resize(nsp_);
  fion.resize(nsp_);
  ftmp.resize(3*namax_);
  
  eself_ = 0.0;
  
  for ( int is = 0; is < nsp_; is++ )
  {
    Species *s = atoms.species_list[is];
    
    const int na = atoms.na(is);
    tau0[is].resize(3*na);
    taum[is].resize(3*na);
    fion_esr[is].resize(3*na);
    fion[is].resize(3*na);
    
    eself_ += na * s->eself();
    na_[is] = na;
    
    zv_[is] = s->zval();
    rcps_[is] = s->rcps();
    
    g = vbasis_->g_ptr();
    double v,dv;  
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhops[is][ig] = s->rhopsg(g[ig]) * omega_inv;
      s->dvlocg(g[ig],v,dv);
      vps[is][ig] =  v * omega_inv;
      dvps[is][ig] =  dv * omega_inv;
    }    
  }
  
  sf.init(tau0,*vbasis_);
  
  cell_moved();
  
  atoms_moved();
  
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::atoms_moved(void)
{
  const AtomSet& atoms = s_.atoms;
  int ngloc = vbasis_->localsize();

  // fill tau0, taum with values in atom_list
  
  atoms.get_positions(tau0);
  sf.update(tau0,*vbasis_);
  
  // compute Fourier coefficients of the local potential
  memset( (void*)&vion_local_g[0], 0, 2*ngloc*sizeof(double) );
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < atoms.nsp(); is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhopst[ig] += s[ig] * rhops[is][ig];
      vion_local_g[ig] += s[ig] * vps[is][ig];
    }
  }
  
  nlp->update_eigr(tau0);
  nlp->update_anl();
  
  // compute esr: pseudocharge repulsion energy
  const UnitCell& cell = s_.wf.cell();
  esr_  = 0.0;
  //desrda[0] = 0.0;
  //desrda[1] = 0.0;
  //desrda[2] = 0.0;
  
  for ( int is = 0; is < nsp_; is++ )
    for ( int i = 0; i < fion_esr[is].size(); i++ )
      fion_esr[is][i] = 0.0;

  for ( int k = 0; k < nsp_; k++ )
  {
    for ( int j = k; j < nsp_; j++ )
    {
      double rckj = sqrt ( rcps_[k]*rcps_[k]+rcps_[j]*rcps_[j] );
      int lax = na_[k];
      if ( k == j ) lax--;

      for ( int l = 0; l < lax; l++ )
      {
        int inf = 0;
        if ( k == j ) inf = l+1;
        for ( int m = inf; m < na_[j]; m++ )
        {
          double xlm = tau0[k][3*l+0] - tau0[j][3*m+0];
          double ylm = tau0[k][3*l+1] - tau0[j][3*m+1];
          double zlm = tau0[k][3*l+2] - tau0[j][3*m+2];
          D3vector vlm(xlm,ylm,zlm);
          cell.fold_in_ws(vlm);
          xlm = vlm.x;
          ylm = vlm.y;
          zlm = vlm.z;
          double rlm = sqrt(xlm*xlm + ylm*ylm + zlm*zlm);
          double arg = rlm / rckj;
          double addesr = zv_[k] * zv_[j] * erfc(arg) / rlm;
          esr_ += addesr;

          double rxlm = vlm.x;
          double rylm = vlm.y;
          double rzlm = vlm.z;
          double addpre = 2.0 * zv_[k]*zv_[j]*exp(-arg*arg)/rckj/sqrt(M_PI);
          double repand = (addesr+addpre) / ( rlm*rlm );
          fion_esr[k][3*l+0] += repand * rxlm;
          fion_esr[j][3*m+0] -= repand * rxlm;
          fion_esr[k][3*l+1] += repand * rylm;
          fion_esr[j][3*m+1] -= repand * rylm;
          fion_esr[k][3*l+2] += repand * rzlm;
          fion_esr[j][3*m+2] -= repand * rzlm;
#if compute_stress
          desrda[0] -= repand * xlm * xlm;
          desrda[1] -= repand * ylm * ylm;
          desrda[2] -= repand * zlm * zlm;
#endif
        }
      }
    }
  }
//   desrda[0] /= al[0];
//   desrda[1] /= al[1];
//   desrda[2] /= al[2];
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::cell_moved(void)
{
  nlp->update_twnl();
}
