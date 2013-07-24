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
// Wavefunction.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Wavefunction.h"
#include "SlaterDet.h"
#include "jacobi.h"
#include "SharedFilePtr.h"
#include <vector>
#include <iomanip>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Context& ctxt) : ctxt_(ctxt), nel_(0),
nempty_(0), nspin_(1), deltaspin_(0), ecut_(0.0), nrowmax_(32)
{
  // create a default wavefunction: one k point, k=0
  kpoint_.resize(1);
  kpoint_[0] = D3vector(0,0,0);
  weight_.resize(1);
  weight_[0] = 1.0;
  compute_nst();
  create_contexts();
  allocate();
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::Wavefunction(const Wavefunction& wf) : ctxt_(wf.ctxt_),
nel_(wf.nel_), nempty_(wf.nempty_), nspin_(wf.nspin_),
deltaspin_(wf.deltaspin_), nrowmax_(wf.nrowmax_),
cell_(wf.cell_), refcell_(wf.refcell_),
ecut_(wf.ecut_), weight_(wf.weight_), kpoint_(wf.kpoint_)
{
  // Create a Wavefunction using the dimensions of the argument

  compute_nst();

  // Next lines: do special allocation of contexts to ensure that
  // contexts are same as those of wf
  spincontext_ = new Context(*wf.spincontext_);
  kpcontext_ = new Context(*wf.kpcontext_);
  sdcontext_ = new Context(*wf.sdcontext_);

  allocate();

  resize(cell_,refcell_,ecut_);
  init();
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::~Wavefunction()
{
  deallocate();
  delete spincontext_;
  delete kpcontext_;
  delete sdcontext_;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::allocate(void)
{
  // create SlaterDets using kpoint list
  const int nkp = kpoint_.size();
  sd_.resize(nspin_);
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    sd_[ispin].resize(nkp);
    for ( int ikp = 0; ikp < nkp; ikp++ )
    {
      sd_[ispin][ikp] = new SlaterDet(*sdcontext_,kpoint_[ikp]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::deallocate(void)
{
  for ( int ispin = 0; ispin < sd_.size(); ispin++ )
  {
    for ( int ikp = 0; ikp < sd_[ispin].size(); ikp++ )
    {
      delete sd_[ispin][ikp];
    }
    sd_[ispin].resize(0);
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
  for ( int ispin = 0; ispin < nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < nkp(); ikp++ )
    {
      sum += weight_[ikp] * sd(ispin,ikp)->entropy(nspin_);
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::resize(const UnitCell& cell, const UnitCell& refcell,
  double ecut)
{
  try
  {
    // resize all SlaterDets using cell, refcell, ecut and nst_[ispin]
    for ( int ispin = 0; ispin < nspin_; ispin++ )
    {
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        sd_[ispin][ikp]->resize(cell,refcell,ecut,nst_[ispin]);
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
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->init();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::clear(void)
{
  // initialize all SlaterDets with zero
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->c().clear();
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
  resize(cell_,refcell_,ecut_);
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
  resize(cell_,refcell_,ecut_);
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
  resize(cell_,refcell_,ecut_);
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
    if ( ctxt_.onpe0() )
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
  resize(cell_,refcell_,ecut_);
  init();
  update_occ(0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::create_contexts(void)
{
  // determine dimensions of sdcontext
  assert(nrowmax_>0);
  const int size = ctxt_.size();
  int npr = nrowmax_;
  while ( size%npr != 0 ) npr--;
  // npr now divides size
  int npc = size/npr;

  spincontext_ = new Context(npr,npc);
  kpcontext_ = new Context(*spincontext_);
  sdcontext_ = new Context(*kpcontext_);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nrowmax(int n)
{
  if ( n > ctxt_.size() )
  {
    if ( ctxt_.onpe0() )
      cout << " Wavefunction::set_nrowmax: nrowmax > ctxt_.size()" << endl;
    return;
  }

  deallocate();
  delete spincontext_;
  delete kpcontext_;
  delete sdcontext_;
  nrowmax_ = n;
  create_contexts();
  compute_nst();
  allocate();
  resize(cell_,refcell_,ecut_);
  init();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::add_kpoint(D3vector kpoint, double weight)
{
  for ( int i = 0; i < kpoint_.size(); i++ )
  {
    if ( length(kpoint - kpoint_[i]) < 1.e-6 )
    {
      if ( ctxt_.onpe0() )
        cout << " Wavefunction::add_kpoint: kpoint already defined"
           << endl;
      return;
    }
  }

  kpoint_.push_back(kpoint);
  weight_.push_back(weight);

  const int nkp = kpoint_.size();
  sd_.resize(nspin_);
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    sd_[ispin].push_back(new SlaterDet(*sdcontext_,kpoint_[nkp-1]));
    sd_[ispin][nkp-1]->resize(cell_,refcell_,ecut_,nst_[ispin]);
  }

  if ( nspin_ == 1 )
  {
    sd_[0][nkp-1]->update_occ(nel_,nspin_);
  }
  else if ( nspin_ == 2 )
  {
    const int nocc_up = (nel_+1)/2+deltaspin_;
    const int nocc_dn = nel_/2 - deltaspin_;
    sd_[0][nkp-1]->update_occ(nocc_up,nspin_);
    sd_[1][nkp-1]->update_occ(nocc_dn,nspin_);
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
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
    if ( ctxt_.onpe0() )
      cout << " Wavefunction::del_kpoint: no such kpoint"
         << endl;
    return;
  }
  deallocate();
  kpoint_.erase(pk);
  weight_.erase(pw);
  allocate();
  resize(cell_,refcell_,ecut_);
  init();
  update_occ(0.0);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize(double amplitude)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->randomize(amplitude);
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
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        sd_[0][ikp]->update_occ(nel_,nspin_);
      }
    }
    else if ( nspin_ == 2 )
    {
      const int nocc_up = (nel_+1)/2+deltaspin_;
      const int nocc_dn = nel_/2 - deltaspin_;
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        sd_[0][ikp]->update_occ(nocc_up,nspin_);
        sd_[1][ikp]->update_occ(nocc_dn,nspin_);
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
    for ( int ispin = 0; ispin < nspin_; ispin++ )
    {
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        sd_[ispin][ikp]->update_occ(nspin_,mu,temp);
        rhosum += weight_[ikp] * sd_[ispin][ikp]->total_charge();
      }
    }

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
      for ( int ispin = 0; ispin < nspin_; ispin++ )
      {
        for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
        {
          sd_[ispin][ikp]->update_occ(nspin_,mu,temp);
          rhosum += weight_[ikp] * sd_[ispin][ikp]->total_charge();
        }
      }
    }

    if ( niter == maxiter )
    {
      cout << "Wavefunction::update_occ: mu did not converge in "
           << maxiter << " iterations" << endl;
      ctxt_.abort(1);
    }

    if ( ctxt_.onpe0() )
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
        for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
        {
          cout << " k = " << kpoint_[ikp] << endl;
          for ( int n = 0; n < sd_[ispin][ikp]->nst(); n++ )
          {
            cout << setw(7) << setprecision(4) << sd_[ispin][ikp]->occ(n);
            if ( ( n%10 ) == 9 ) cout << endl;
          }
          if ( sd_[ispin][ikp]->nst() % 10 != 0 )
            cout << endl;
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::gram(void)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->gram();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::riccati(Wavefunction& wf)
{
  assert(wf.context() == ctxt_);
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->riccati(*wf.sd_[ispin][ikp]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::align(Wavefunction& wf)
{
  assert(wf.context() == ctxt_);
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sd_[ispin][ikp]->align(*wf.sd_[ispin][ikp]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
complex<double> Wavefunction::dot(const Wavefunction& wf) const
{
  assert(wf.context() == ctxt_);
  complex<double> sum = 0.0;
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      sum += sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::diag(Wavefunction& dwf, bool eigvec)
{
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
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print(ostream& os, string encoding, string tag) const
{
  if ( ctxt_.onpe0() )
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

  if ( ctxt_.onpe0() )
    os << "</" << tag << ">" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write(SharedFilePtr& sfp, string encoding, string tag) const
{
  sfp.sync();

  if ( ctxt_.onpe0() )
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

  if ( ctxt_.onpe0() )
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
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::info(ostream& os, string tag) const
{
  if ( ctxt_.onpe0() )
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

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( ctxt_.onpe0() )
        cout << " kpoint: " << kpoint_[ikp] << " weight: " << weight_[ikp]
             << endl;
      sd_[ispin][ikp]->info(os);
    }
  }

  if ( ctxt_.onpe0() )
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
  assert(ctxt_ == wf.ctxt_);
  assert(nel_ == wf.nel_);
  assert(nempty_== wf.nempty_);
  assert(nspin_ == wf.nspin_);
  assert(nrowmax_ == wf.nrowmax_);
  assert(deltaspin_ == wf.deltaspin_);
  //!!assert(cell_ == wf.cell_);
  assert(refcell_ == wf.refcell_);
  assert(ecut_ == wf.ecut_);

  cell_ = wf.cell_;
  weight_ = wf.weight_;
  kpoint_ = wf.kpoint_;

  assert(*spincontext_ == *wf.spincontext_);
  assert(*kpcontext_   == *wf.kpcontext_);
  assert(*sdcontext_   == *wf.sdcontext_);

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      *sd_[ispin][ikp] = *wf.sd_[ispin][ikp];
    }
  }
  return *this;
}
