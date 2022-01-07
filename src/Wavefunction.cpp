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
#include "cout0.h"
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
  //cout << MPIdata::rank() << ": nkp=" << nkp << endl;
  //cout << MPIdata::rank() << ": MPIdata::nkpb()=" << MPIdata::nkpb() << endl;
  for ( int ikpb = 0; ikpb < MPIdata::nkpb(); ++ikpb )
  {
    nkp_loc_[ikpb] = nkp / MPIdata::nkpb() +
                     (ikpb < (nkp % MPIdata::nkpb()) ? 1 : 0);
    //cout << MPIdata::rank() << ": nkp_loc_[" << ikpb << "]="
    //     << nkp_loc_[ikpb] << endl;
  }

  // round robin allocation of kpoints
  ikp_global_.resize(0);
  for ( int ikpg = MPIdata::ikpb(); ikpg < nkp; ikpg += MPIdata::nkpb() )
  {
    ikp_global_.push_back(ikpg);
  }

  //cout << MPIdata::rank() << ": nkp_loc_[MPIdata::ikpb]="
  //     << nkp_loc_[MPIdata::ikpb()] << endl;
  //cout << MPIdata::rank() << ": ikp_global_.size()="
  //     << ikp_global_.size() << endl;
  assert(nkp_loc_[MPIdata::ikpb()] == ikp_global_.size());

  // compute local number of spins nsp_loc_[ispb]
  nsp_loc_.resize(MPIdata::nspb());
  for ( int ispb = 0; ispb < MPIdata::nspb(); ++ispb)
  {
    nsp_loc_[ispb] = nspin_ / MPIdata::nspb() +
                      (ispb < (nspin_ % MPIdata::nspb()) ? 1 : 0);
    //cout << MPIdata::rank() << ": nsp_loc_[" << ispb << "]="
    //     << nsp_loc_[ispb] << endl;
  }

  // round robin allocation of spins
  isp_global_.resize(0);
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
int Wavefunction::nkp_loc() const { return nkp_loc_[MPIdata::ikpb()]; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::ikp_global(int ikp) const { return ikp_global_[ikp]; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::ikp_local(int ikpg) const
{
  // return local index of ikpg if hosted on this task, -1 otherwise
  if ( ( ikpg % MPIdata::nkpb() ) == MPIdata::ikpb() )
    return ikpg / MPIdata::nkpb();
  else
    return -1;
}

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
int Wavefunction::nsp_loc() const { return nsp_loc_[MPIdata::ispb()]; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::isp_global(int isp) const { return isp_global_[isp]; }

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::isp_local(int ispg) const
{
  // return local index of ispg if hosted on this task, -1 otherwise
  if ( ( ispg % MPIdata::nspb() ) == MPIdata::ispb() )
    return ispg / MPIdata::nspb();
  else
    return -1;
}

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
  double tsum;
  MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::kp_sp_comm());
  return tsum;
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

  // Add SlaterDet
  // determine on which block the SlaterDet should be added
  const int ikp_new = kpoint_.size() - 1;
  sd_.resize(nsp_loc());
  const int ikp_loc = ikp_local(ikp_new);
  if ( ikp_loc >= 0 )
  {
    // new kpoint is local to the current block
    nkp_loc_[MPIdata::ikpb()]++;
    ikp_global_.push_back(ikp_new);
    for ( int isp_loc = 0; isp_loc < nsp_loc(); ++isp_loc )
    {
      sd_[isp_loc].push_back(new SlaterDet(sd_ctxt_,kpoint_[ikp_new]));
      sd_[isp_loc][ikp_loc]->resize(cell_,refcell_,ecut_,
        nst_[isp_global(isp_loc)]);
    }

    if ( nspin_ == 1 )
    {
      if ( nsp_loc() > 0 )
        sd_[0][ikp_loc]->update_occ(nel_,nspin_);
    }
    else if ( nspin_ == 2 )
    {
      for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
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
    else
    {
      // incorrect value of nspin_
      assert(false);
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
  // check if new kpoint position coincides with existing kpoint
  for ( int i = 0; i < kpoint_.size(); i++ )
  {
    if ( length(newkpoint - kpoint_[i]) < 1.e-6 )
    {
      if ( MPIdata::onpe0() )
        cout << " Wavefunction::move_kpoint: kpoint already defined "
             << "at kpoint new position"
             << endl;
      return;
    }
  }

  // ikp: global index of kpoint to be moved
  const int ikp_loc = ikp_local(ikp);

  // copy wavefunctions from old SlaterDet at kpoint to new SlaterDet
  // at newkpoint
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    const int isp_loc = isp_local(ispin);
    if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
    {
      // this task holds kpoint ikp

      // create new SlaterDet at newkpoint
      SlaterDet *sd = sd_[isp_loc][ikp_loc];
      SlaterDet *sdn = new SlaterDet(sd->context(),newkpoint);
      sdn->resize(cell_,refcell_,ecut_,nst_[ispin]);
      sdn->init();
      // copy wave functions from old to new SlaterDet
      const Basis& basis = sd_[isp_loc][ikp_loc]->basis();
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
      delete sd_[isp_loc][ikp_loc];
      // reassign pointer
      sd_[isp_loc][ikp_loc] = sdn;

      // randomize new SlaterDet to break symmetry
      sd_[isp_loc][ikp_loc]->randomize(1.e-4);

    } // if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
  } // for ispin

  kpoint_[ikp] = newkpoint;

  // update occupation numbers
  if ( nspin_ == 1 )
  {
    if ( ( isp_local(0) >= 0 ) && ( ikp_loc >= 0 ) )
    {
      sd_[isp_local(0)][ikp_loc]->update_occ(nel_,nspin_);
    }
  }
  else if ( nspin_ == 2 )
  {
    const int nocc_up = (nel_+1)/2+deltaspin_;
    const int nocc_dn = nel_/2 - deltaspin_;

    if ( ( isp_local(0) >= 0 ) && ( ikp_loc >= 0 ) )
      sd_[isp_local(0)][ikp_loc]->update_occ(nocc_up,nspin_);

    if ( ( isp_local(1) >= 0 ) && ( ikp_loc >= 0 ) )
      sd_[isp_local(1)][ikp_loc]->update_occ(nocc_dn,nspin_);
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize(double amplitude)
{
  srand48((long int) MPIdata::rank());
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
    MPI_Allreduce(&rhosum,&tmpsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::kp_sp_comm());
    rhosum = tmpsum;

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
      MPI_Allreduce(&rhosum,&tmpsum,1,MPI_DOUBLE,MPI_SUM,MPIdata::kp_sp_comm());
      rhosum = tmpsum;
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
      cout << "<mu> " << setprecision(6) << mu / eVolt << " </mu>" << endl;
      cout << "<occ_set>" << endl;
    }

    for ( int ispin = 0; ispin < nspin(); ++ispin )
    {
      const int isp_loc = isp_local(ispin);
      for ( int ikp = 0; ikp < nkp(); ++ikp )
      {
        const int ikp_loc = ikp_local(ikp);
        ostringstream ostr;
        ostr.setf(ios::fixed,ios::floatfield);
        ostr.setf(ios::right,ios::adjustfield);
        int isrc = -1;
        if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
        {
          if ( MPIdata::sd_rank() == 0 )
          {
            ostr.str("");
            isrc = MPIdata::rank();
            const int nst = sd(isp_loc,ikp_loc)->nst();
            ostr <<    "  <occ spin=\"" << ispin
                 << "\" kpoint=\""
                 << setprecision(8)
                 << sd(isp_loc,ikp_loc)->kpoint()
                 << "\" weight=\""
                 << setprecision(8)
                 << weight(ikp)
                 << "\" n=\"" << nst << "\">" << endl;
            for ( int i = 0; i < nst; i++ )
            {
              ostr << setw(7) << setprecision(4)
                   << sd(isp_loc,ikp_loc)->occ(i);
              if ( i%10 == 9 ) ostr << endl;
            }
            if ( nst%10 != 0 ) ostr << endl;
            ostr << "  </occ>" << endl;
          }
        }
        cout0(ostr.str(),isrc);
        MPI_Barrier(MPIdata::comm());
      }
    }
    if ( MPIdata::onpe0() )
       cout << "</occ_set>" << endl;
    MPI_Barrier(MPIdata::comm());
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
  for ( int isp_loc = 0; isp_loc < sd_.size(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < sd_[isp_loc].size(); ++ikp_loc )
    {
      const int ikpg = ikp_global(ikp_loc);
      sum += weight_[ikpg] *
        sd_[isp_loc][ikp_loc]->dot(*wf.sd_[isp_loc][ikp_loc]);
    }
  }
  // sum over kpoint and spin comms
  complex<double> tsum;
  MPI_Allreduce(&sum,&tsum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPIdata::kp_sp_comm());
  return tsum;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::diag(Wavefunction& dwf, bool eigvec)
{
  // subspace diagonalization of <*this | dwf>
  // if eigvec==true, eigenvectors are computed and stored in *this, dwf is
  // overwritten
  for ( int isp_loc = 0; isp_loc < nsp_loc(); ++isp_loc )
  {
    assert( nst_[isp_global(isp_loc)] > 0 );
    for ( int ikp_loc = 0; ikp_loc < nkp_loc(); ++ikp_loc )
    {
      // compute eigenvalues
      if ( sd(isp_loc,ikp_loc)->basis().real() )
      {
        // proxy real matrices c, cp
        DoubleMatrix c(sd(isp_loc,ikp_loc)->c());
        DoubleMatrix cp(dwf.sd(isp_loc,ikp_loc)->c());

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
        sd(isp_loc,ikp_loc)->set_eig(w);
      }
      else
      {
        ComplexMatrix& c(sd(isp_loc,ikp_loc)->c());
        ComplexMatrix& cp(dwf.sd(isp_loc,ikp_loc)->c());
        ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
        h.gemm('c','n',1.0,c,cp,0.0);
        // cout << " Hamiltonian at k = "
        //      << sd(ispin,ikp)->kpoint() << endl;
        // cout << h;
        valarray<double> w(h.m());
        if ( eigvec )
        {
          ComplexMatrix z(c.context(),c.n(),c.n(),c.nb(),c.nb());
          h.heev('l',w,z);
          cp = c;
          c.gemm('n','n',1.0,cp,z,0.0);
        }
        else
        {
          h.heev('l',w);
        }
        // set eigenvalues in SlaterDet
        sd(isp_loc,ikp_loc)->set_eig(w);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print(ostream& os, string encoding, string tag) const
{
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

  for ( int ispin = 0; ispin < nspin(); ++ispin )
  {
    const int isp_loc = isp_local(ispin);
    for ( int ikp = 0; ikp < nkp(); ++ikp )
    {
      const int ikp_loc = ikp_local(ikp);
      MPI_Barrier(MPIdata::comm());
      if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
      {
        sd_[isp_loc][ikp_loc]->print(os,encoding,weight_[ikp],ispin,nspin_);
      }
    }
  }
  MPI_Barrier(MPIdata::comm());

  if ( MPIdata::onpe0() )
    os << "</" << tag << ">" << endl;
}

///////////////////////////////////////////////////////////////////////////////
void Wavefunction::write(SharedFilePtr& sfp, string encoding, string tag) const
{
  assert(sizeof(size_t)==sizeof(MPI_Offset));
  assert(sizeof(long long int)==sizeof(MPI_Offset));

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
    size_t len = str.size();
    MPI_Status status;
    int err = MPI_File_write_at(sfp.file(),sfp.offset(),(void*)str.c_str(),
              len,MPI_CHAR,&status);
    if ( err != 0 )
      cout << " Wavefunction::write: error in MPI_File_write" << endl;
    sfp.advance(len);
  }

  sfp.sync();

  vector<vector<string> > sdstr;
  sdstr.resize(nsp_loc());
  for ( int isp_loc = 0; isp_loc < nsp_loc(); ++isp_loc )
  {
    const int ispin = isp_global(isp_loc);
    sdstr[isp_loc].resize(nkp_loc());
    for ( int ikp_loc = 0; ikp_loc < nkp_loc(); ++ikp_loc )
    {
      const int ikp = ikp_global(ikp_loc);
      // serialize sd[isp_loc][ikp_loc] into sdstr[isp_loc][ikp_loc]
      sd_[isp_loc][ikp_loc]->str(sdstr[isp_loc][ikp_loc],encoding,
                                 weight_[ikp],ispin,nspin_);
    }
  }

  // sdstr now contains all data to be written
  // write data in order of increasing ispin, ikp
  for ( int ispin = 0; ispin < nspin(); ++ispin )
  {
    for ( int ikp = 0; ikp < nkp(); ++ikp )
    {
      MPI_Barrier(MPIdata::comm());
      const int isp_loc = isp_local(ispin);
      const int ikp_loc = ikp_local(ikp);
      const char* wbuf = 0;
      size_t len = 0;
      MPI_Offset off = 0;
      if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
      {
        wbuf = sdstr[isp_loc][ikp_loc].c_str();
        len = sdstr[isp_loc][ikp_loc].size();
        // compute offset of data on current task
        long long int local_offset = 0;
        long long int local_size = len;
        MPI_Scan(&local_size, &local_offset, 1,
                 MPI_LONG_LONG, MPI_SUM, MPIdata::sd_comm());
        // correct for inclusive scan by subtracting local_size
        local_offset -= local_size;
        off = sfp.offset() + local_offset;
      }
      MPI_Status status;
      int err = MPI_File_write_at_all(sfp.file(),off,wbuf,len,
                                      MPI_CHAR,&status);
      int count;
      MPI_Get_count(&status,MPI_CHAR,&count);
      assert(count==len);

      if ( err != 0 )
        cout << " Wavefunction::write: error in MPI_File_write_at_all" << endl;
      sfp.set_offset(off+len);
      sfp.sync();
    }
  }

  if ( MPIdata::onpe0() )
  {
    ostringstream os;
    os << "</" << tag << ">" << endl;
    string str(os.str());
    size_t len = str.size();
    MPI_Status status;
    int err = MPI_File_write_at(sfp.file(),sfp.offset(),(void*)str.c_str(),
              len,MPI_CHAR,&status);
    if ( err != 0 )
      cout << " Wavefunction::write: error in MPI_File_write" << endl;
    sfp.advance(len);
  }
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

  for ( int ispin = 0; ispin < nspin(); ++ispin )
  {
    const int isp_loc = isp_local(ispin);
    for ( int ikp = 0; ikp < nkp(); ++ikp )
    {
      const int ikp_loc = ikp_local(ikp);
      string s;
      int isrc = -1;
      if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
      {
        s = sd_[isp_loc][ikp_loc]->info();
        if ( MPIdata::sd_rank() == 0 )
          isrc = MPIdata::rank();
      }
      cout0(s,isrc);
      MPI_Barrier(MPIdata::comm());
    }
  }

  if ( MPIdata::onpe0() )
    os << "</" << tag << ">" << endl;
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
