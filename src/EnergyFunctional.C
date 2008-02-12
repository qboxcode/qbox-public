////////////////////////////////////////////////////////////////////////////////
//
// EnergyFunctional.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: EnergyFunctional.C,v 1.29 2008-02-12 05:39:19 fgygi Exp $

#include "EnergyFunctional.h"
#include "Sample.h"
#include "Species.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "StructureFactor.h"
#include "XCPotential.h"
#include "NonLocalPotential.h"
#include "ConfinementPotential.h"

#include "Timer.h"
#include "blas.h"

#include <iostream>
#include <iomanip>
#include <algorithm> // fill()
using namespace std;

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::EnergyFunctional(const Sample& s, const ChargeDensity& cd)
 : s_(s), cd_(cd)
{
  const AtomSet& atoms = s_.atoms;
  const Wavefunction& wf = s_.wf;

  sigma_ekin.resize(6);
  sigma_econf.resize(6);
  sigma_eps.resize(6);
  sigma_ehart.resize(6);
  sigma_exc.resize(6);
  sigma_enl.resize(6);
  sigma_esr.resize(6);
  sigma.resize(6);

  vbasis_ = cd_.vbasis();
  //cout << vbasis_->context().mype() << ": vbasis_->context() = "
  //     << vbasis_->context() << endl;

  // define FT's on vbasis contexts

  vft = cd_.vft();
  int np0v = vft->np0();
  int np1v = vft->np1();
  int np2v = vft->np2();

  v_r.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    v_r[ispin].resize(vft->np012loc());
  }
  tmp_r.resize(vft->np012loc());

  if ( s_.ctxt_.onpe0() )
  {
    cout << "  EnergyFunctional: np0v,np1v,np2v: " << np0v << " "
         << np1v << " " << np2v << endl;
    cout << "  EnergyFunctional: vft->np012(): "
         << vft->np012() << endl;
  }

  const int ngloc = vbasis_->localsize();
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
  nlp.resize(wf.nkp());
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    nlp[ikp] = new NonLocalPotential(s_.atoms, *wf.sd(0,ikp));
  }

  vion_local_g.resize(ngloc);
  dvion_local_g.resize(ngloc);
  vlocal_g.resize(ngloc);
  vtemp.resize(ngloc);
  rhoelg.resize(ngloc);
  rhogt.resize(ngloc);
  rhopst.resize(ngloc);

  tau0.resize(nsp_);
  fion_esr.resize(nsp_);
  ftmp.resize(3*namax_);

  eself_ = 0.0;

  for ( int is = 0; is < nsp_; is++ )
  {
    Species *s = atoms.species_list[is];

    const int na = atoms.na(is);
    tau0[is].resize(3*na);
    fion_esr[is].resize(3*na);

    eself_ += na * s->eself();
    na_[is] = na;

    zv_[is] = s->zval();
    rcps_[is] = s->rcps();
  }

  // FT for interpolation of wavefunctions on the fine grid
  ft.resize(wf.nkp());
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    ft[ikp] = cd_.ft(ikp);
  }

  // Confinement potentials
  cfp.resize(wf.nkp());
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    cfp[ikp] = 0;
    const double facs = 2.0;
    const double sigmas = 0.5;
    cfp[ikp] =
      new ConfinementPotential(s_.ctrl.ecuts,facs,sigmas,
        wf.sd(0,ikp)->basis());
  }

  sf.init(tau0,*vbasis_);

  cell_moved();

  atoms_moved();

}

////////////////////////////////////////////////////////////////////////////////
EnergyFunctional::~EnergyFunctional(void)
{
  delete xcp;
  for ( int ikp = 0; ikp < s_.wf.nkp(); ikp++ )
  {
    delete cfp[ikp];
    delete nlp[ikp];
  }

  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<timing name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::update_vhxc(void)
{
  // called when the charge density has changed
  // update Hartree and xc potentials using the charge density cd_
  // compute Hartree and xc energies

  const Wavefunction& wf = s_.wf;
  const UnitCell& cell = wf.cell();
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const double *const g2i = vbasis_->g2i_ptr();
  const double fpi = 4.0 * M_PI;
  const int ngloc = vbasis_->localsize();
  double tsum[2];

  // compute total electronic density: rhoelg = rho_up + rho_dn
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

  // update XC energy and potential
  tmap["exc"].start();
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    fill(v_r[ispin].begin(),v_r[ispin].end(),0.0);
  xcp->update(v_r);
  exc_ = xcp->exc();
  tmap["exc"].stop();

  // compute local potential energy:
  // integral of el. charge times ionic local pot.
  int len=2*ngloc,inc1=1;
  tsum[0] = 2.0 * ddot(&len,(double*)&rhoelg[0],&inc1,
         (double*)&vion_local_g[0],&inc1);
  // remove double counting for G=0
  if ( vbasis_->context().myrow() == 0 )
  {
    tsum[0] -= real(conj(rhoelg[0])*vion_local_g[0]);
  }
  tsum[0] *= omega; // tsum[0] contains eps

  // Hartree energy
  ehart_ = 0.0;
  double ehsum = 0.0;
  for ( int ig = 0; ig < ngloc; ig++ )
  {
    const complex<double> tmp = rhoelg[ig] + rhopst[ig];
    ehsum += norm(tmp) * g2i[ig];
    rhogt[ig] = tmp;
  }
  // factor 1/2 from definition of Ehart cancels with half sum over G
  // Note: rhogt[ig] includes a factor 1/Omega
  // Factor omega in next line yields prefactor 4 pi / omega in
  tsum[1] = omega * fpi * ehsum;
  // tsum[1] contains ehart

  vbasis_->context().dsum(2,1,&tsum[0],2);
  eps_   = tsum[0];
  ehart_ = tsum[1];

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
}

////////////////////////////////////////////////////////////////////////////////
double EnergyFunctional::energy(bool compute_hpsi, Wavefunction& dwf,
              bool compute_forces, vector<vector<double> >& fion,
              bool compute_stress, valarray<double>& sigma)
{
  const bool debug_stress = compute_stress &&
    s_.ctrl.debug.find("STRESS") != string::npos;
  const double fpi = 4.0 * M_PI;

  const Wavefunction& wf = s_.wf;
  const UnitCell& cell(wf.cell());
  const double omega = cell.volume();
  const double omega_inv = 1.0 / omega;
  const int ngloc = vbasis_->localsize();
  const double *const g2i = vbasis_->g2i_ptr();
  assert(wf.nspin()==1);

  const bool use_confinement = s_.ctrl.ecuts > 0.0;

  for ( int is = 0; is < fion.size(); is++ )
    for ( int ia = 0; ia < fion[is].size(); ia++ )
      fion[is][ia] = 0.0;

  if ( compute_hpsi )
  {
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
        dwf.sd(ispin,ikp)->c().clear();
  }

  // kinetic energy
  tmap["ekin"].start();

  // compute ekin, confinement energy, stress from ekin and econf
  // ekin = sum_G |c_G|^2  G^2
  // econf = sum_G |c_G|^2 fstress[G]
  // stress_ekin_ij = (1/Omega) sum_G |c_G|^2 * 2 * G_i * G_j
  // stress_econf_ij = (1/Omega) sum_G |c_G|^2 * dfstress[G] * G_i * G_j
  ekin_ = 0.0;
  econf_ = 0.0;
  sigma_ekin = 0.0;
  sigma_econf = 0.0;
  valarray<double> sum(0.0,14), tsum(0.0,14);
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      const double weight = wf.weight(ikp);
      const SlaterDet& sd = *(wf.sd(ispin,ikp));
      const Basis& wfbasis = sd.basis();
      const D3vector kp = wfbasis.kpoint();
      // factor fac in next lines: 2.0 for G and -G (if basis is real) and
      // 0.5 from 1/(2m)
      const double fac = wfbasis.real() ? 1.0 : 0.5;
      const ComplexMatrix& c = sd.c();
      const Context& sdctxt = sd.context();

      // compute psi2sum(G) = fac * sum_G occ(n) psi2(n,G)
      const int ngwloc = wfbasis.localsize();
      valarray<double> psi2sum(ngwloc);
      const complex<double>* p = c.cvalptr();
      const int mloc = c.mloc();
      const int nloc = c.nloc();
      // nn = global n index
      const int nnbase = sdctxt.mycol() * c.nb();
      const double * const occ = sd.occ_ptr(nnbase);
      for ( int ig = 0; ig < ngwloc; ig++ )
      {
        double tmpsum = 0.0;
        for ( int n = 0; n < nloc; n++ )
        {
          const double psi2 = norm(p[ig+n*mloc]);
          tmpsum += fac * occ[n] * psi2;
        }
        psi2sum[ig] = tmpsum;
      }

      // accumulate contributions to ekin,econf,sigma_ekin,sigma_econf in tsum
      const double *const kpg2  = wfbasis.kpg2_ptr();
      const double *const kpg_x = wfbasis.kpgx_ptr(0);
      const double *const kpg_y = wfbasis.kpgx_ptr(1);
      const double *const kpg_z = wfbasis.kpgx_ptr(2);
      tsum = 0.0;

      for ( int ig = 0; ig < ngwloc; ig++ )
      {
        const double psi2s = psi2sum[ig];

        // tsum[0]: ekin partial sum
        tsum[0] += psi2s * kpg2[ig];

        if ( compute_stress )
        {
          const double tkpgx = kpg_x[ig];
          const double tkpgy = kpg_y[ig];
          const double tkpgz = kpg_z[ig];

          const double fac_ekin = 2.0 * psi2s;

          tsum[1]  += fac_ekin * tkpgx * tkpgx;
          tsum[2]  += fac_ekin * tkpgy * tkpgy;
          tsum[3]  += fac_ekin * tkpgz * tkpgz;
          tsum[4]  += fac_ekin * tkpgx * tkpgy;
          tsum[5]  += fac_ekin * tkpgy * tkpgz;
          tsum[6]  += fac_ekin * tkpgx * tkpgz;

        }
        // tsum[0-6] contains the contributions to
        // ekin, sigma_ekin, from vector ig
      } // ig

      if ( use_confinement )
      {
        const valarray<double>& fstress = cfp[ikp]->fstress();
        const valarray<double>& dfstress = cfp[ikp]->dfstress();
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          const double psi2s = psi2sum[ig];
          // tsum[7]: econf partial sum
          tsum[7] += psi2s * fstress[ig];

          if ( compute_stress )
          {
            const double tkpgx = kpg_x[ig];
            const double tkpgy = kpg_y[ig];
            const double tkpgz = kpg_z[ig];

            const double fac_econf = psi2s * dfstress[ig];
            tsum[8]  += fac_econf * tkpgx * tkpgx;
            tsum[9]  += fac_econf * tkpgy * tkpgy;
            tsum[10] += fac_econf * tkpgz * tkpgz;
            tsum[11] += fac_econf * tkpgx * tkpgy;
            tsum[12] += fac_econf * tkpgy * tkpgz;
            tsum[13] += fac_econf * tkpgx * tkpgz;
          }
          // tsum[7-13] contains the contributions to
          // econf,sigma_econf from vector ig
        } // ig
      }

      sum += weight * tsum;
    } // ikp
  } // ispin

  // sum contains the contributions to ekin, etc.. from this task
  wf.context().dsum(14,1,&sum[0],14);

  ekin_  = sum[0];
  sigma_ekin[0] = sum[1];
  sigma_ekin[1] = sum[2];
  sigma_ekin[2] = sum[3];
  sigma_ekin[3] = sum[4];
  sigma_ekin[4] = sum[5];
  sigma_ekin[5] = sum[6];

  econf_ = sum[7];
  sigma_econf[0] = sum[8];
  sigma_econf[1] = sum[9];
  sigma_econf[2] = sum[10];
  sigma_econf[3] = sum[11];
  sigma_econf[4] = sum[12];
  sigma_econf[5] = sum[13];

  sigma_ekin *= omega_inv;
  sigma_econf *= omega_inv;

  tmap["ekin"].stop();

  // Stress from Eps
  sigma_eps = 0.0;
  if ( compute_stress )
  {
    tsum = 0.0;
    const double *const gi  = vbasis_->gi_ptr();
    const double *const g_x = vbasis_->gx_ptr(0);
    const double *const g_y = vbasis_->gx_ptr(1);
    const double *const g_z = vbasis_->gx_ptr(2);
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      // factor of 2 in next line: G and -G
      // note: gi[0] == 0.0 in next line (no division by zero)
      const double fac = 2.0 * gi[ig] *
        real( conj(rhoelg[ig]) * dvion_local_g[ig] );

      const double tgx = g_x[ig];
      const double tgy = g_y[ig];
      const double tgz = g_z[ig];

      tsum[0] += fac * tgx * tgx;
      tsum[1] += fac * tgy * tgy;
      tsum[2] += fac * tgz * tgz;
      tsum[3] += fac * tgx * tgy;
      tsum[4] += fac * tgy * tgz;
      tsum[5] += fac * tgx * tgz;
    }
    vbasis_->context().dsum(6,1,&tsum[0],6);

    sigma_eps[0] = eps_ * omega_inv + tsum[0];
    sigma_eps[1] = eps_ * omega_inv + tsum[1];
    sigma_eps[2] = eps_ * omega_inv + tsum[2];
    sigma_eps[3] = tsum[3];
    sigma_eps[4] = tsum[4];
    sigma_eps[5] = tsum[5];
  }

  // Stress from Hartree energy
  if ( compute_stress )
  {
    tsum = 0.0;
    const double *const g_x = vbasis_->gx_ptr(0);
    const double *const g_y = vbasis_->gx_ptr(1);
    const double *const g_z = vbasis_->gx_ptr(2);

    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const double temp = norm(rhogt[ig]) * g2i[ig] * g2i[ig];
      const double tgx = g_x[ig];
      const double tgy = g_y[ig];
      const double tgz = g_z[ig];

      tsum[0] += temp * tgx * tgx;
      tsum[1] += temp * tgy * tgy;
      tsum[2] += temp * tgz * tgz;
      tsum[3] += temp * tgx * tgy;
      tsum[4] += temp * tgy * tgz;
      tsum[5] += temp * tgx * tgz;
    }

    for ( int is = 0; is < nsp_; is++ )
    {
      // Factor 2 in next line: 2*real(x) = ( x + h.c. )
      // Factor 0.5 in next line: definition of sigma
      const double fac2 = 2.0 * 0.25 * rcps_[is]*rcps_[is];
      complex<double> *s = &sf.sfac[is][0];
      for ( int ig = 0; ig < ngloc; ig++ )
      {
        const complex<double> sg = s[ig];
        const complex<double> rg = rhogt[ig];
        // next line: keep only real part
        const double temp = fac2 *
        ( rg.real() * sg.real() +
          rg.imag() * sg.imag() )
        * rhops[is][ig] * g2i[ig];

        const double tgx = g_x[ig];
        const double tgy = g_y[ig];
        const double tgz = g_z[ig];

        tsum[0] += temp * tgx * tgx;
        tsum[1] += temp * tgy * tgy;
        tsum[2] += temp * tgz * tgz;
        tsum[3] += temp * tgx * tgy;
        tsum[4] += temp * tgy * tgz;
        tsum[5] += temp * tgx * tgz;
      }
    }
    vbasis_->context().dsum(6,1,&tsum[0],6);
    // Factor in next line:
    //  factor 2 from G and -G
    //  factor fpi from definition of sigma
    //  no factor 1/Omega^2 (already included in rhogt[ig] above)
    sigma_ehart[0] = ehart_ * omega_inv - 2.0 * fpi * tsum[0];
    sigma_ehart[1] = ehart_ * omega_inv - 2.0 * fpi * tsum[1];
    sigma_ehart[2] = ehart_ * omega_inv - 2.0 * fpi * tsum[2];
    sigma_ehart[3] = - 2.0 * fpi * tsum[3];
    sigma_ehart[4] = - 2.0 * fpi * tsum[4];
    sigma_ehart[5] = - 2.0 * fpi * tsum[5];
  } // compute_stress

  // Stress from exchange-correlation
  if ( compute_stress )
  {
    xcp->compute_stress(sigma_exc);
  }

  // Non local energy and forces
  tmap["nonlocal"].start();
  // modify next loop to span only local ikp
  enl_ = 0.0;
  vector<vector<double> > fion_enl;
  fion_enl.resize(nsp_);
  for ( int is = 0; is < nsp_; is++ )
    fion_enl[is].resize(3*na_[is]);
  valarray<double> sigma_enl_kp(6);
  sigma_enl = 0.0;
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    enl_ += wf.weight(ikp) * nlp[ikp]->energy(compute_hpsi,*dwf.sd(0,ikp),
            compute_forces, fion_enl, compute_stress, sigma_enl_kp);

    if ( compute_forces )
      for ( int is = 0; is < nsp_; is++ )
        for ( int i = 0; i < 3*na_[is]; i++ )
          fion[is][i] += wf.weight(ikp) * fion_enl[is][i];

    if ( compute_stress )
      sigma_enl += wf.weight(ikp) * sigma_enl_kp;
  }
  // add here sum of enl across rows of kpcontext
  tmap["nonlocal"].stop();

  ecoul_ = ehart_ + esr_ - eself_;
  ets_ = 0.0;
  if ( s_.ctrl.fermi_temp > 0.0 )
  {
    const double wf_entropy = wf.entropy();
    const double boltz = 1.0 / ( 11605.0 * 2.0 * 13.6058 );
    ets_ = - wf_entropy * s_.ctrl.fermi_temp * boltz;
  }
  etotal_ = ekin_ + econf_ + eps_ + enl_ + ecoul_ + exc_ + ets_;

  if ( compute_hpsi )
  {
    tmap["hpsi"].start();
    assert(wf.nspin()==1);
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        const SlaterDet& sd = *(wf.sd(ispin,ikp));
        SlaterDet& sdp = *(dwf.sd(ispin,ikp));
        const ComplexMatrix& c = sd.c();
        const Basis& wfbasis = sd.basis();
        ComplexMatrix& cp = dwf.sd(ispin,ikp)->c();
        const int mloc = cp.mloc();
        const double* kpg2 = wfbasis.kpg2_ptr();
        const int ngwloc = wfbasis.localsize();

        // Laplacian
        if ( use_confinement )
        {
          for ( int n = 0; n < sd.nstloc(); n++ )
          {
            const valarray<double>& fstress = cfp[ikp]->fstress();
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              cp[ig+mloc*n] += 0.5 * ( kpg2[ig] + fstress[ig] ) *
                               c[ig+mloc*n];
            }
          }
        }
        else
        {
          for ( int n = 0; n < sd.nstloc(); n++ )
          {
            for ( int ig = 0; ig < ngwloc; ig++ )
            {
              cp[ig+mloc*n] += 0.5 * kpg2[ig] * c[ig+mloc*n];
            }
          }
        }
        sd.rs_mul_add(*ft[ikp], &v_r[ispin][0], sdp);
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
          // complex<double> teigr = ei0[kv[3*ig+0]] *
          //                         ei1[kv[3*ig+1]] *
          //                         ei2[kv[3*ig+2]];
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

  sigma = sigma_ekin + sigma_econf + sigma_eps + sigma_enl +
          sigma_ehart + sigma_exc + sigma_esr;
  if ( debug_stress && s_.ctxt_.onpe0() )
  {
    const double gpa = 29421.5;
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout << setprecision(8);
    cout << " <stress_tensor unit=\"atomic_units\">\n"
         << "   <sigma_ekin_xx> " << setw(12)
         << sigma_ekin[0] << " </sigma_ekin_xx>\n"
         << "   <sigma_ekin_yy> " << setw(12)
         << sigma_ekin[1] << " </sigma_ekin_yy>\n"
         << "   <sigma_ekin_zz> " << setw(12)
         << sigma_ekin[2] << " </sigma_ekin_zz>\n"
         << "   <sigma_ekin_xy> " << setw(12)
         << sigma_ekin[3] << " </sigma_ekin_xy>\n"
         << "   <sigma_ekin_yz> " << setw(12)
         << sigma_ekin[4] << " </sigma_ekin_yz>\n"
         << "   <sigma_ekin_xz> " << setw(12)
         << sigma_ekin[5] << " </sigma_ekin_xz>\n"
         << endl
         << "   <sigma_econf_xx> " << setw(12)
         << sigma_econf[0] << " </sigma_econf_xx>\n"
         << "   <sigma_econf_yy> " << setw(12)
         << sigma_econf[1] << " </sigma_econf_yy>\n"
         << "   <sigma_econf_zz> " << setw(12)
         << sigma_econf[2] << " </sigma_econf_zz>\n"
         << "   <sigma_econf_xy> " << setw(12)
         << sigma_econf[3] << " </sigma_econf_xy>\n"
         << "   <sigma_econf_yz> " << setw(12)
         << sigma_econf[4] << " </sigma_econf_yz>\n"
         << "   <sigma_econf_xz> " << setw(12)
         << sigma_econf[5] << " </sigma_econf_xz>\n"
         << endl
         << "   <sigma_eps_xx> " << setw(12)
         << sigma_eps[0] << " </sigma_eps_xx>\n"
         << "   <sigma_eps_yy> " << setw(12)
         << sigma_eps[1] << " </sigma_eps_yy>\n"
         << "   <sigma_eps_zz> " << setw(12)
         << sigma_eps[2] << " </sigma_eps_zz>\n"
         << "   <sigma_eps_xy> " << setw(12)
         << sigma_eps[3] << " </sigma_eps_xy>\n"
         << "   <sigma_eps_yz> " << setw(12)
         << sigma_eps[4] << " </sigma_eps_yz>\n"
         << "   <sigma_eps_xz> " << setw(12)
         << sigma_eps[5] << " </sigma_eps_xz>\n"
         << endl
         << "   <sigma_enl_xx> " << setw(12)
         << sigma_enl[0] << " </sigma_enl_xx>\n"
         << "   <sigma_enl_yy> " << setw(12)
         << sigma_enl[1] << " </sigma_enl_yy>\n"
         << "   <sigma_enl_zz> " << setw(12)
         << sigma_enl[2] << " </sigma_enl_zz>\n"
         << "   <sigma_enl_xy> " << setw(12)
         << sigma_enl[3] << " </sigma_enl_xy>\n"
         << "   <sigma_enl_yz> " << setw(12)
         << sigma_enl[4] << " </sigma_enl_yz>\n"
         << "   <sigma_enl_xz> " << setw(12)
         << sigma_enl[5] << " </sigma_enl_xz>\n"
         << endl
         << "   <sigma_ehart_xx> " << setw(12)
         << sigma_ehart[0] << " </sigma_ehart_xx>\n"
         << "   <sigma_ehart_yy> " << setw(12)
         << sigma_ehart[1] << " </sigma_ehart_yy>\n"
         << "   <sigma_ehart_zz> " << setw(12)
         << sigma_ehart[2] << " </sigma_ehart_zz>\n"
         << "   <sigma_ehart_xy> " << setw(12)
         << sigma_ehart[3] << " </sigma_ehart_xy>\n"
         << "   <sigma_ehart_yz> " << setw(12)
         << sigma_ehart[4] << " </sigma_ehart_yz>\n"
         << "   <sigma_ehart_xz> " << setw(12)
         << sigma_ehart[5] << " </sigma_ehart_xz>\n"
         << endl
         << "   <sigma_exc_xx> " << setw(12)
         << sigma_exc[0] << " </sigma_exc_xx>\n"
         << "   <sigma_exc_yy> " << setw(12)
         << sigma_exc[1] << " </sigma_exc_yy>\n"
         << "   <sigma_exc_zz> " << setw(12)
         << sigma_exc[2] << " </sigma_exc_zz>\n"
         << "   <sigma_exc_xy> " << setw(12)
         << sigma_exc[3] << " </sigma_exc_xy>\n"
         << "   <sigma_exc_yz> " << setw(12)
         << sigma_exc[4] << " </sigma_exc_yz>\n"
         << "   <sigma_exc_xz> " << setw(12)
         << sigma_exc[5] << " </sigma_exc_xz>\n"
         << endl
         << "   <sigma_esr_xx> " << setw(12)
         << sigma_esr[0] << " </sigma_esr_xx>\n"
         << "   <sigma_esr_yy> " << setw(12)
         << sigma_esr[1] << " </sigma_esr_yy>\n"
         << "   <sigma_esr_zz> " << setw(12)
         << sigma_esr[2] << " </sigma_esr_zz>\n"
         << "   <sigma_esr_xy> " << setw(12)
         << sigma_esr[3] << " </sigma_esr_xy>\n"
         << "   <sigma_esr_yz> " << setw(12)
         << sigma_esr[4] << " </sigma_esr_yz>\n"
         << "   <sigma_esr_xz> " << setw(12)
         << sigma_esr[5] << " </sigma_esr_xz>\n"
         << endl
         << "   <sigma_eks_xx> " << setw(12) << sigma[0] << " </sigma_eks_xx>\n"
         << "   <sigma_eks_yy> " << setw(12) << sigma[1] << " </sigma_eks_yy>\n"
         << "   <sigma_eks_zz> " << setw(12) << sigma[2] << " </sigma_eks_zz>\n"
         << "   <sigma_eks_xy> " << setw(12) << sigma[3] << " </sigma_eks_xy>\n"
         << "   <sigma_eks_yz> " << setw(12) << sigma[4] << " </sigma_eks_yz>\n"
         << "   <sigma_eks_xz> " << setw(12) << sigma[5] << " </sigma_eks_xz>\n"
         << " </stress_tensor>" << endl;
  }


  return etotal_;
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::atoms_moved(void)
{
  const AtomSet& atoms = s_.atoms;
  int ngloc = vbasis_->localsize();

  // fill tau0 with values in atom_list

  atoms.get_positions(tau0);
  sf.update(tau0,*vbasis_);

  // compute Fourier coefficients of the local potential
  memset( (void*)&vion_local_g[0], 0, 2*ngloc*sizeof(double) );
  memset( (void*)&dvion_local_g[0], 0, 2*ngloc*sizeof(double) );
  memset( (void*)&rhopst[0], 0, 2*ngloc*sizeof(double) );
  for ( int is = 0; is < atoms.nsp(); is++ )
  {
    complex<double> *s = &sf.sfac[is][0];
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      const complex<double> sg = s[ig];
      rhopst[ig] += sg * rhops[is][ig];
      vion_local_g[ig] += sg * vps[is][ig];
      dvion_local_g[ig] += sg * dvps[is][ig];
    }
  }

  // compute esr: pseudocharge repulsion energy
  const UnitCell& cell = s_.wf.cell();
  const double omega_inv = 1.0 / cell.volume();

  esr_  = 0.0;
  sigma_esr = 0.0;
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
          double esr_lm = zv_[k] * zv_[j] * erfc(arg) / rlm;
          esr_ += esr_lm;

          double desr_erfc = 2.0 * zv_[k]*zv_[j]*exp(-arg*arg)/rckj/sqrt(M_PI);

          // desrdr = (1/r) d Esr / dr
          double desrdr = -(esr_lm+desr_erfc) / ( rlm*rlm );
          fion_esr[k][3*l+0] -= desrdr * xlm;
          fion_esr[j][3*m+0] += desrdr * xlm;
          fion_esr[k][3*l+1] -= desrdr * ylm;
          fion_esr[j][3*m+1] += desrdr * ylm;
          fion_esr[k][3*l+2] -= desrdr * zlm;
          fion_esr[j][3*m+2] += desrdr * zlm;

          sigma_esr[0] += desrdr * xlm * xlm;
          sigma_esr[1] += desrdr * ylm * ylm;
          sigma_esr[2] += desrdr * zlm * zlm;
          sigma_esr[3] += desrdr * xlm * ylm;
          sigma_esr[4] += desrdr * ylm * zlm;
          sigma_esr[5] += desrdr * xlm * zlm;
        }
      }
    }
  }
  sigma_esr *= - omega_inv;
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::cell_moved(void)
{
  const Wavefunction& wf = s_.wf;
  const UnitCell& cell = wf.cell();
  // resize vbasis_
  vbasis_->resize(cell,s_.wf.refcell(),4.0*s_.wf.ecut());

  const int ngloc = vbasis_->localsize();
  const double omega = cell.volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  const AtomSet& atoms = s_.atoms;
  for ( int is = 0; is < nsp_; is++ )
  {
    Species *s = atoms.species_list[is];
    const double * const g = vbasis_->g_ptr();
    double v,dv;
    for ( int ig = 0; ig < ngloc; ig++ )
    {
      rhops[is][ig] = s->rhopsg(g[ig]) * omega_inv;
      s->dvlocg(g[ig],v,dv);
      vps[is][ig] =  v * omega_inv;
      dvps[is][ig] =  dv * omega_inv;
    }
  }

  // Update confinement potentials and non-local potentials
  for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
  {
    cfp[ikp]->update();
    nlp[ikp]->update_twnl();
  }
}

////////////////////////////////////////////////////////////////////////////////
void EnergyFunctional::print(ostream& os) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setprecision(8);
  os << "  <ekin>   " << setw(15) << ekin() << " </ekin>\n"
     << "  <econf>  " << setw(15) << econf() << " </econf>\n"
     << "  <eps>    " << setw(15) << eps() << " </eps>\n"
     << "  <enl>    " << setw(15) << enl() << " </enl>\n"
     << "  <ecoul>  " << setw(15) << ecoul() << " </ecoul>\n"
     << "  <exc>    " << setw(15) << exc() << " </exc>\n"
     << "  <esr>    " << setw(15) << esr() << " </esr>\n"
     << "  <eself>  " << setw(15) << eself() << " </eself>\n"
     << "  <ets>    " << setw(15) << ets() << " </ets>\n"
     << "  <etotal> " << setw(15) << etotal() << " </etotal>\n"
     << flush;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const EnergyFunctional& e )
{
  e.print(os);
  return os;
}
