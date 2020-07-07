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
// Wavefunction.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "MPIdata.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "jacobi.h"
#include "SharedFilePtr.h"
#include <vector>
#include <iomanip>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Context& sd_ctxt) : sd_ctxt_(sd_ctxt), nel_(0),
nempty_(0), nspin_(1), deltaspin_(0), ecut_(0.0)
{
  // create a default wavefunction: one k point, k=0
  kpoint_.resize(1);
  kpoint_[0] = D3vector(0,0,0);
  weight_.resize(1);
  weight_[0] = 1.0;
  compute_nst();
  allocate();
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Wavefunction& wf) : sd_ctxt_(wf.sd_ctxt_),
nel_(wf.nel_), nempty_(wf.nempty_), nspin_(wf.nspin_), nst_(wf.nst_),
deltaspin_(wf.deltaspin_), cell_(wf.cell_), refcell_(wf.refcell_),
ecut_(wf.ecut_), weight_(wf.weight_), kpoint_(wf.kpoint_)
{
  allocate();
  resize();
  init();
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::~Wavefunction()
{
  deallocate();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::allocate(void)
{
  // compute local number of kpoints nkp_loc_[ikpb]
  nkp_loc_.resize(MPIdata::nkpb());
  const int nkp = kpoint_.size();
  cout << MPIdata::rank() << ": nkp=" << nkp << endl;
  cout << MPIdata::rank() << ": MPIdata::nkpb()=" << MPIdata::nkpb() << endl;
  for ( int ikpb = 0; ikpb < MPIdata::nkpb(); ++ikpb )
  {
    nkp_loc_[ikpb] = nkp / MPIdata::nkpb() +
                     (ikpb < (nkp % MPIdata::nkpb()) ? 1 : 0);
    cout << MPIdata::rank() << ": nkp_loc_[" << ikpb << "]="
         << nkp_loc_[ikpb] << endl;
  }

  // round robin allocation of kpoints
  for ( int ikpg = MPIdata::ikpb(); ikpg < nkp; ikpg += MPIdata::nkpb() )
  {
    ikp_global_.push_back(ikpg);
    cout << MPIdata::rank() << ": ikp_global_.push_back(" << ikpg << ")\n";
  }

  cout << MPIdata::rank() << ": nkp_loc_[MPIdata::ikpb]="
       << nkp_loc_[MPIdata::ikpb()] << endl;
  cout << MPIdata::rank() << ": ikp_global_.size()="
       << ikp_global_.size() << endl;
  assert(nkp_loc_[MPIdata::ikpb()] == ikp_global_.size());

  // compute local number of spins nsp_loc_[ispb]
  nsp_loc_.resize(MPIdata::nspb());
  for ( int ispb = 0; ispb < MPIdata::nspb(); ++ispb)
  {
    nsp_loc_[ispb] = nspin_ / MPIdata::nspb() +
                      (ispb < (nspin_ % MPIdata::nspb()) ? 1 : 0);
    cout << MPIdata::rank() << ": nsp_loc_[" << ispb << "]="
         << nsp_loc_[ispb] << endl;
  }

  // round robin allocation of spins
  for ( int ispg = MPIdata::ispb(); ispg < nspin_; ispg += MPIdata::nspb() )
  {
    isp_global_.push_back(ispg);
  }

  assert(nsp_loc_[MPIdata::ispb()] == isp_global_.size());

  // create local SlaterDets
  sd_.resize(nsp_loc_[MPIdata::ispb()]);
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    sd_[isp_loc].resize(nkp_loc_[MPIdata::ikpb()]);
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      int ikp = ikp_global_[ikp_loc];
      sd_[isp_loc][ikp_loc] = new SlaterDet(sd_ctxt_,kpoint_[ikp]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::deallocate(void)
{
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      delete sd_[isp_loc][ikp_loc];
    }
    sd_[isp_loc].resize(0);
  }
}

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nkp(void) const { return kpoint_.size(); }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nel() const { return nel_; } // total number of electrons

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nst() const
{
  if ( nspin_ == 1 )
    return nst_[0];
  else
    return nst_[0]+nst_[1];
}

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nst(int ispin) const
{
  assert(ispin >= 0 && ispin < 2);
  return nst_[ispin];
}

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nempty() const { return nempty_; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nspin() const { return nspin_; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::deltaspin() const { return deltaspin_; }

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::entropy(void) const
{
  double sum = 0.0;
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      int ikp = ikp_global_[ikp_loc];
      sum += weight_[ikp] * sd(isp_loc,ikp_loc)->entropy(nspin_);
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::resize(void)
{
  try
  {
    // resize all SlaterDets using current cell_, refcell_, ecut_, nst_[ispin]
    for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
      {
        int isp = isp_global_[isp_loc];
        sd_[isp_loc][ikp_loc]->resize(cell_,refcell_,ecut_,nst_[isp]);
      }
    }
  }
  catch ( const SlaterDetException& sdex )
  {
    cout << " Wavefunction: SlaterDetException during resize: " << endl
         << sdex.msg << endl;
    // no resize took place
    return;
  }
  catch ( bad_alloc )
  {
    cout << " Wavefunction: insufficient memory for resize operation" << endl;
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::resize(const UnitCell& cell, const UnitCell& refcell,
  double ecut)
{
  try
  {
    // resize all SlaterDets using cell, refcell, ecut and nst_[ispin]
    for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
      {
        int isp = isp_global_[isp_loc];
        sd_[isp_loc][ikp_loc]->resize(cell,refcell,ecut,nst_[isp]);
      }
    }
    cell_ = cell;
    refcell_ = refcell;
    ecut_ = ecut;
  }
  catch ( const SlaterDetException& sdex )
  {
    cout << " Wavefunction: SlaterDetException during resize: " << endl
         << sdex.msg << endl;
    // no resize took place
    return;
  }
  catch ( bad_alloc )
  {
    cout << " Wavefunction: insufficient memory for resize operation" << endl;
    return;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::init(void)
{
  // initialize all SlaterDets with lowest energy plane waves
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->init();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::clear(void)
{
  // initialize all SlaterDets with zero
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->c().clear();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::reset(void)
{
  // reset to single kpoint, ecut=0
  deallocate();
  nel_ = 0;
  nempty_ = 0;
  nspin_ = 1;
  deltaspin_ = 0;
  kpoint_.resize(1);
  kpoint_[0] = D3vector(0,0,0);
  weight_.resize(1);
  weight_[0] = 1.0;
  compute_nst();
  allocate();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::compute_nst(void)
{
  // recompute nst from nel_, deltaspin_, nempty_

  nst_.resize(nspin_);
  if ( nspin_ == 1 )
  {
    nst_[0] = ( nel_ + 1 ) / 2 + nempty_;
  }
  else
  {
    // nspin == 2
    nst_[0] = ( nel_ + 1 ) / 2 + deltaspin_ + nempty_;
    nst_[1] = nel_ / 2 - deltaspin_ + nempty_;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nel(int nel)
{
  if ( nel == nel_ ) return;
  if ( nel < 0 )
  {
    cout << " Wavefunction::set_nel: nel < 0" << endl;
    return;
  }

  nel_ = nel;
  compute_nst();
  resize();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nempty(int nempty)
{
  if ( nempty == nempty_ ) return;
  if ( nempty < 0 )
  {
    cout << " Wavefunction::set_nempty: negative value" << endl;
    return;
  }
  nempty_ = nempty;
  compute_nst();
  update_occ(0.0);
  resize();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nspin(int nspin)
{
  assert(nspin==1 || nspin==2);
  if ( nspin == nspin_ ) return;

  deallocate();
  nspin_ = nspin;

  compute_nst();
  allocate();
  resize();
  init();
  update_occ(0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_deltaspin(int deltaspin)
{
  if ( deltaspin == deltaspin_ ) return;

  // check if value of deltaspin would result in nst[1] < 0
  // nst_[1] = nel_ / 2 - deltaspin_ + nempty_;
  if ( nel_ / 2 - deltaspin + nempty_ < 0 )
  {
    if ( MPIdata::onpe0() )
    {
      cout << " Wavefunction::set_deltaspin: nel+nempty too small" << endl;
      cout << " Wavefunction::set_deltaspin: cannot set deltaspin" << endl;
      return;
    }
  }
  deallocate();
  deltaspin_ = deltaspin;
  compute_nst();
  allocate();
  resize();
  init();
  update_occ(0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::add_kpoint(D3vector kpoint, double weight)
{
  for ( int i = 0; i < kpoint_.size(); i++ )
  {
    if ( length(kpoint - kpoint_[i]) < 1.e-6 )
    {
      if ( MPIdata::onpe0() )
        cout << " Wavefunction::add_kpoint: kpoint already defined"
             << endl;
      return;
    }
  }

  kpoint_.push_back(kpoint);
  weight_.push_back(weight);

  // determine on which kpoint block the next kpoint should be added
  const int nkp = kpoint_.size();
  const int ikpbnew = (nkp+1) % MPIdata::nkpb();
  if ( ikpbnew == MPIdata::ikpb() )
  {
    sd_.resize(nspin_);
    for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
    {
      sd_[isp_loc].push_back(new SlaterDet(sd_ctxt_,kpoint_[nkp-1]));
      int ispin = isp_global_[isp_loc];
      sd_[isp_loc][nkp-1]->resize(cell_,refcell_,ecut_,nst_[ispin]);
      if ( nspin_ == 1 )
      {
        sd_[0][nkp-1]->update_occ(nel_,nspin_);
      }
      else if ( nspin_ == 2 )
      {
        if ( isp_global_[isp_loc] == 0 )
        {
          const int nocc_up = (nel_+1)/2+deltaspin_;
          sd_[isp_loc][nkp-1]->update_occ(nocc_up,nspin_);
        }
        else
        {
          const int nocc_dn = nel_/2 - deltaspin_;
          sd_[isp_loc][nkp-1]->update_occ(nocc_dn,nspin_);
        }
      }
      else
      {
        // incorrect value of nspin_
        assert(false);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::del_kpoint(D3vector kpoint)
{
  bool found = false;
  vector<D3vector>::iterator pk = kpoint_.begin();
  vector<double>::iterator pw = weight_.begin();
  while ( !found && pk != kpoint_.end() )
  {
    if ( length(kpoint - *pk) < 1.e-6 )
    {
      found = true;
    }
    else
    {
      pk++;
      pw++;
    }
  }
  if ( !found )
  {
    if ( MPIdata::onpe0() )
      cout << " Wavefunction::del_kpoint: no such kpoint"
         << endl;
    return;
  }
  deallocate();
  kpoint_.erase(pk);
  weight_.erase(pw);
  allocate();
  resize();
  init();
  update_occ(0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::move_kpoint(D3vector kpoint, D3vector newkpoint)
{
  bool found = false;
  int ikp = 0;
  while ( !found && ikp < kpoint_.size() )
  {
    if ( length(kpoint_[ikp] - kpoint) < 1.e-6 )
    {
      found = true;
    }
    else
    {
      ikp++;
    }
  }
  if ( !found )
  {
    if ( MPIdata::onpe0() )
      cout << " Wavefunction::move_kpoint: no such kpoint"
         << endl;
    return;
  }
  // check if newkpoint coincides with existing kpoint
  for ( int i = 0; i < kpoint_.size(); i++ )
  {
    if ( length(newkpoint - kpoint_[i]) < 1.e-6 )
    {
      if ( MPIdata::onpe0() )
        cout << " Wavefunction::move_kpoint: kpoint already defined "
             << "at newkpoint position"
             << endl;
      return;
    }
  }

  cout << "Wavefunction::move_kpoint: not implemented" << endl;
#if 0
  // copy wavefunctions from old SlaterDet at kpoint to new SlaterDet
  // at newkpoint
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    // create new SlaterDet at newkpoint
    SlaterDet *sd = sd_[ispin][ikp];
    SlaterDet *sdn = new SlaterDet(*sdcontext_,newkpoint);
    sdn->resize(cell_,refcell_,ecut_,nst_[ispin]);
    sdn->init();
    // copy wave functions from old to new SlaterDet
    const Basis& basis = sd_[ispin][ikp]->basis();
    const Basis& newbasis = sdn->basis();
    // transform all states to real space and interpolate
    int np0 = max(basis.np(0),newbasis.np(0));
    int np1 = max(basis.np(1),newbasis.np(1));
    int np2 = max(basis.np(2),newbasis.np(2));
    FourierTransform ft1(basis,np0,np1,np2);
    FourierTransform ft2(newbasis,np0,np1,np2);
    // allocate real-space grid
    valarray<complex<double> > tmpr(ft1.np012loc());
    for ( int n = 0; n < sd->nstloc(); n++ )
    {
      for ( int i = 0; i < tmpr.size(); i++ )
        tmpr[i] = 0.0;
      ComplexMatrix& c = sd->c();
      ComplexMatrix& cn = sdn->c();
      ft1.backward(c.valptr(n*c.mloc()),&tmpr[0]);
      // if the new kpoint is Gamma, remove the phase of the wavefunction
      // using the phase of the first element of the array
      if ( newbasis.real() )
      {
        double arg0 = arg(tmpr[0]);
        // broadcast argument found on task 0
        MPI_Bcast(&arg0,1,MPI_DOUBLE,0,basis.comm());
        complex<double> z = polar(1.0,-arg0);
        for ( int i = 0; i < tmpr.size(); i++ )
          tmpr[i] *= z;
      }
      ft2.forward(&tmpr[0], cn.valptr(n*cn.mloc()));
    }

    // delete old SlaterDet
    delete sd_[ispin][ikp];
    // reassign pointer
    sd_[ispin][ikp] = sdn;
  }
  kpoint_[ikp] = newkpoint;

  if ( nspin_ == 1 )
  {
    sd_[0][ikp]->update_occ(nel_,nspin_);
  }
  else if ( nspin_ == 2 )
  {
    const int nocc_up = (nel_+1)/2+deltaspin_;
    const int nocc_dn = nel_/2 - deltaspin_;
    sd_[0][ikp]->update_occ(nocc_up,nspin_);
    sd_[1][ikp]->update_occ(nocc_dn,nspin_);
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize(double amplitude)
{
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->randomize(amplitude);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::update_occ(double temp)
{
  // update occupation numbers using eigenvalues in SlaterDet
  if ( temp == 0.0 )
  {
    // zero temperature
    if ( nspin_ == 1 )
    {
      for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
        {
          sd_[isp_loc][ikp_loc]->update_occ(nel_,nspin_);
        }
      }
    }
    else if ( nspin_ == 2 )
    {
      for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
        {
          if ( isp_global_[isp_loc] == 0 )
          {
            const int nocc_up = (nel_+1)/2+deltaspin_;
            sd_[isp_loc][ikp_loc]->update_occ(nocc_up,nspin_);
          }
          else
          {
            const int nocc_dn = nel_/2 - deltaspin_;
            sd_[isp_loc][ikp_loc]->update_occ(nocc_dn,nspin_);
          }
        }
      }
    }
    else
    {
      // incorrect value of nspin_
      assert(false);
    }
  }
  else
  {
    // finite temperature
    const double eVolt = 0.036749023; // 1 eV in Hartree
    const int maxiter = 500;

    // loop to find value of mu
    double mu = 0.0;
    double dmu = 2.0 * eVolt;
    const double totalcharge = (double) nel_;
    enum direction { up, down };
    direction dir = up;

    double rhosum = 0.0;
    for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
    {
      for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
      {
        sd_[isp_loc][ikp_loc]->update_occ(nspin_,mu,temp);
        const int ikpg = ikp_global_[ikp_loc];
        rhosum += weight_[ikpg] * sd_[isp_loc][ikp_loc]->total_charge();
      }
    }
    // sum over spin and kpoints
    double tmpsum = 0.0;
    MPI_Allreduce(&rhosum,&tmpsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::kp_comm());
    MPI_Allreduce(&tmpsum,&rhosum,1,MPI_DOUBLE,MPI_SUM,MPIdata::sp_comm());

    int niter = 0;
    while ( niter < maxiter && fabs(rhosum - nel_) > 1.e-10 )
    {
      niter++;
      if ( rhosum < totalcharge )
      {
        if ( dir == down ) dmu /= 2.0;
        mu += dmu;
        dir = up;
      }
      else
      {
        if ( dir == up ) dmu /= 2.0;
        mu -= dmu;
        dir = down;
      }
      rhosum = 0.0;
      for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
      {
        for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
        {
          sd_[isp_loc][ikp_loc]->update_occ(nspin_,mu,temp);
          int ikp = ikp_global_[ikp_loc];
          rhosum += weight_[ikp] * sd_[isp_loc][ikp_loc]->total_charge();
        }
      }
      MPI_Allreduce(&rhosum,&tmpsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::kp_comm());
      MPI_Allreduce(&tmpsum,&rhosum,1,MPI_DOUBLE,MPI_SUM,MPIdata::sp_comm());
    }

    if ( niter == maxiter )
    {
      cout << "Wavefunction::update_occ: mu did not converge in "
           << maxiter << " iterations" << endl;
      MPI_Abort(MPIdata::comm(),1);
    }

    if ( MPIdata::onpe0() )
    {
      cout << " Wavefunction::update_occ: sum = "
           << rhosum << endl;
      cout << " Wavefunction::update_occ: mu = "
           << setprecision(4) << mu / eVolt << " eV" << endl;

      cout.setf(ios::right,ios::adjustfield);
      cout.setf(ios::fixed,ios::floatfield);

      cout << " Wavefunction::update_occ: occupation numbers" << endl;
      for ( int ispin = 0; ispin < nspin_; ispin++ )
      {
        const int isp_loc = ispin / MPIdata::nspb();
        const int ispb = ispin % MPIdata::nspb();
        if ( ispb == MPIdata::ispb() )
        {
          for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
          {
            const int ikp_loc = ikp / MPIdata::nkpb();
            const int ikpb = ikp % MPIdata::nkpb();
            if ( ikpb == MPIdata::ikpb() )
            {
              cout << " k = " << kpoint_[ikp] << endl;
              for ( int n = 0; n < sd_[isp_loc][ikp_loc]->nst(); n++ )
              {
                cout << setw(7) << setprecision(4)
                     << sd_[isp_loc][ikp_loc]->occ(n);
                if ( ( n%10 ) == 9 ) cout << endl;
              }
              if ( sd_[ispin][ikp]->nst() % 10 != 0 )
                cout << endl;
            }
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::gram(void)
{
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->gram();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::riccati(Wavefunction& wf)
{
  assert(wf.sd_context() == sd_ctxt_);
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->riccati(*wf.sd_[isp_loc][ikp_loc]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::align(Wavefunction& wf)
{
  assert(wf.sd_context() == sd_ctxt_);
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      sd_[isp_loc][ikp_loc]->align(*wf.sd_[isp_loc][ikp_loc]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
complex<double> Wavefunction::dot(const Wavefunction& wf) const
{
  assert(wf.sd_context() == sd_ctxt_);
  complex<double> sum = 0.0;
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sum += weight_[ikp] * sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::diag(Wavefunction& dwf, bool eigvec)
{
  cout << "Wavefunction::diag: not implemented" << endl;
#if 0
  // subspace diagonalization of <*this | dwf>
  // if eigvec==true, eigenvectors are computed and stored in *this, dwf is
  // overwritten
  for ( int ispin = 0; ispin < nspin(); ispin++ )
  {
    if ( nst_[ispin] > 0 )
    {
      for ( int ikp = 0; ikp < nkp(); ikp++ )
      {
        // compute eigenvalues
        if ( sd(ispin,ikp)->basis().real() )
        {
          // proxy real matrices c, cp
          DoubleMatrix c(sd(ispin,ikp)->c());
          DoubleMatrix cp(dwf.sd(ispin,ikp)->c());

          DoubleMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());

          // factor 2.0 in next line: G and -G
          h.gemm('t','n',2.0,c,cp,0.0);
          // rank-1 update correction
          h.ger(-1.0,c,0,cp,0);

          // cout << " Hamiltonian at k = " << sd(ispin,ikp)->kpoint()
          //      << endl;
          // cout << h;

#if 1
          valarray<double> w(h.m());
          if ( eigvec )
          {
            DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
            h.syev('l',w,z);
            //h.syevx('l',w,z,1.e-6);
            cp = c;
            c.gemm('n','n',1.0,cp,z,0.0);
          }
          else
          {
            h.syev('l',w);
          }
#else
          vector<double> w(h.m());
          DoubleMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
          const int maxsweep = 30;
          int nsweep = jacobi(maxsweep,1.e-6,h,z,w);
          if ( eigvec )
          {
            cp = c;
            c.gemm('n','n',1.0,cp,z,0.0);
          }
#endif
          // set eigenvalues in SlaterDet
          sd(ispin,ikp)->set_eig(w);
        }
        else
        {
          ComplexMatrix& c(sd(ispin,ikp)->c());
          ComplexMatrix& cp(dwf.sd(ispin,ikp)->c());
          ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
          h.gemm('c','n',1.0,c,cp,0.0);
          // cout << " Hamiltonian at k = "
          //      << sd(ispin,ikp)->kpoint() << endl;
          // cout << h;
          valarray<double> w(h.m());
          if ( eigvec )
          {
#if DEBUG
            ComplexMatrix hcopy(h);
#endif
            ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
            h.heev('l',w,z);
            cp = c;
            c.gemm('n','n',1.0,cp,z,0.0);
#if DEBUG
            // check that z contains eigenvectors of h
            // diagonal matrix with eigenvalues on diagonal
            ComplexMatrix d(c.context(),c.n(),c.n(),c.nb(),c.nb());
            // the following test works only on one task
            assert(ctxt_.size()==1);
            for ( int i = 0; i < d.m(); i++ )
              d[i+d.n()*i] = w[i];
            ComplexMatrix dz(c.context(),c.n(),c.n(),c.nb(),c.nb());
            dz.gemm('n','c',1.0,d,z,0.0);
            ComplexMatrix zdz(c.context(),c.n(),c.n(),c.nb(),c.nb());
            zdz.gemm('n','n',1.0,z,dz,0.0);
            // zdz should be equal to hcopy
            zdz -= hcopy;
            cout << " heev: norm of error: " << zdz.nrm2() << endl;
#endif
          }
          else
          {
            h.heev('l',w);
          }
          // set eigenvalues in SlaterDet
          sd(ispin,ikp)->set_eig(w);
        }
      }
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print(ostream& os, string encoding, string tag) const
{
  cout << "Wavefunction::print: not implemented" << endl;
#if 0
  if ( MPIdata::onpe0() )
  {
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os << setprecision(10);
    os << "<domain a=\""
       << cell_.a(0) << "\"\n        b=\""
       << cell_.a(1) << "\"\n        c=\""
       << cell_.a(2) << "\"/>" << endl;
    os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
       <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
       <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
  }

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      sd_[ispin][ikp]->print(os,encoding,weight_[ikp],ispin,nspin_);
  }

  if ( MPIdata::onpe0() )
    os << "</" << tag << ">" << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write(SharedFilePtr& sfp, string encoding, string tag) const
{
  cout << "Wavefunction::print: not implemented" << endl;
#if 0
  sfp.sync();

  if ( MPIdata::onpe0() )
  {
    ostringstream os;
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os << setprecision(10);
    os << "<domain a=\""
       << cell_.a(0) << "\"\n        b=\""
       << cell_.a(1) << "\"\n        c=\""
       << cell_.a(2) << "\"/>" << endl;
    if ( refcell_.volume() != 0.0 )
    {
      os << "<reference_domain a=\""
         << refcell_.a(0) << "\"\n        b=\""
         << refcell_.a(1) << "\"\n        c=\""
         << refcell_.a(2) << "\"/>" << endl;
    }
    os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
       <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
       <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
    string str(os.str());
    int len = str.size();
#if USE_MPI
    MPI_Status status;
    int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)str.c_str(),
              len,MPI_CHAR,&status);
    if ( err != 0 )
      cout << " Wavefunction::write: error in MPI_File_write" << endl;
    sfp.advance(len);
#else
    sfp.file().write(str.c_str(),len);
#endif
  }

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->write(sfp,encoding,weight_[ikp],ispin,nspin_);
    }
  }

  sfp.sync();

  if ( MPIdata::onpe0() )
  {
    ostringstream os;
    os << "</" << tag << ">" << endl;
    string str(os.str());
    int len = str.size();
#if USE_MPI
    MPI_Status status;
    int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*)str.c_str(),
              len,MPI_CHAR,&status);
    if ( err != 0 )
      cout << " Wavefunction::write: error in MPI_File_write" << endl;
    sfp.advance(len);
#else
    sfp.file().write(str.c_str(),len);
#endif
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::info(ostream& os, string tag) const
{
  if ( MPIdata::onpe0() )
  {
    os << "<" << tag << " ecut=\"" << ecut_ << "\""
       << " nspin=\"" << nspin_ << "\""
       << " nel=\"" << nel_ << "\""
       << " nempty=\"" << nempty_ << "\">" << endl;
    os.setf(ios::fixed,ios::floatfield);
    os << "<cell a=\""
       << setprecision(6) << cell_.a(0) << "\"\n      b=\""
       << cell_.a(1) << "\"\n      c=\""
       << cell_.a(2) << "\"/>" << endl;
    os << " reciprocal lattice vectors" << endl
       << setprecision(6)
       << " " << cell_.b(0) << endl
       << " " << cell_.b(1) << endl
       << " " << cell_.b(2) << endl;
    os << "<refcell a=\""
       << refcell_.a(0) << "\"\n         b=\""
       << refcell_.a(1) << "\"\n         c=\""
       << refcell_.a(2) << "\"/>" << endl;
    os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
       <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
       <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
  }

  cout << "Wavefunction::info: not fully implemented" << endl;
#if 0
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      if ( MPIdata::onpe0() )
        cout << " kpoint: " << kpoint_[ikp] << " weight: " << weight_[ikp]
             << endl;
      sd_[ispin][ikp]->info(os);
    }
  }

  if ( MPIdata::onpe0() )
    os << "</" << tag << ">" << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, const Wavefunction& wf)
{
  wf.print(os,"text","wavefunction");
  return os;
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction& Wavefunction::operator=(const Wavefunction& wf)
{
  if ( this == &wf ) return *this;
  assert(sd_ctxt_ == wf.sd_ctxt_);
  assert(nel_ == wf.nel_);
  assert(nempty_== wf.nempty_);
  assert(nspin_ == wf.nspin_);
  assert(deltaspin_ == wf.deltaspin_);
  assert(refcell_ == wf.refcell_);
  assert(ecut_ == wf.ecut_);

  cell_ = wf.cell_;
  weight_ = wf.weight_;
  kpoint_ = wf.kpoint_;

  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      *sd_[isp_loc][ikp_loc] = *wf.sd_[isp_loc][ikp_loc];
    }
  }
  return *this;
}
