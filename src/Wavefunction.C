////////////////////////////////////////////////////////////////////////////////
//
// Wavefunction.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Wavefunction.C,v 1.20 2005-02-04 21:58:59 fgygi Exp $

#include "Wavefunction.h"
#include "SlaterDet.h"
#include <vector>
#include <iomanip>
#if USE_CSTDIO_LFS
#include <sstream>
#include <cstdio>
#endif
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
  
  // create sd contexts and SlaterDets
  const int nkp = kpoint_.size();
  
  assert(ctxt_.active());
  
  assert(nspin_ == 1);
  assert(nkp == 1);
  
  spincontext_.resize(1);
  sdcontext_.resize(1);
  sdcontext_[0].resize(1);
  sd_.resize(1);
  sd_[0].resize(1);
  
  spincontext_[0] = new Context(*wf.spincontext_[0]);
  sdcontext_[0][0] = 0;
  sd_[0][0] = 0;
  if ( spincontext_[0]->active() )
  {
    sdcontext_[0][0] = new Context(*wf.spincontext_[0]);
    sd_[0][0] = new SlaterDet(*sdcontext_[0][0],kpoint_[0]);
  }
  
  resize(cell_,refcell_,ecut_);
  reset();
}

////////////////////////////////////////////////////////////////////////////////
Wavefunction::~Wavefunction()
{
  deallocate();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::allocate(void)
{
  // create sd contexts and SlaterDets
  const int nkp = kpoint_.size();
  
  assert(ctxt_.active());
  
  assert(nspin_ == 1);
  assert(nkp == 1);
  
  // determine dimensions of sdcontext
  assert(nrowmax_>0);
  const int size = ctxt_.size();
  int npr = nrowmax_;
  while ( size%npr != 0 ) npr--;
  // npr now divides size
  int npc = size/npr;

  spincontext_.resize(1);
  sdcontext_.resize(1);
  sdcontext_[0].resize(1);
  sd_.resize(1);
  sd_[0].resize(1);
  
  spincontext_[0] = new Context(npr,npc);
  //cout << *spincontext_[0];
  sdcontext_[0][0] = 0;
  sd_[0][0] = 0;
  if ( spincontext_[0]->active() )
  {
    sdcontext_[0][0] = new Context(*spincontext_[0]);
    sd_[0][0] = new SlaterDet(*sdcontext_[0][0],kpoint_[0]);
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::deallocate(void)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    if ( spincontext_[ispin] != 0 )
    {
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        if ( sdcontext_[ispin][ikp] != 0 )
        {
          delete sd_[ispin][ikp];
          delete sdcontext_[ispin][ikp];
        }
      }
      delete spincontext_[ispin];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::clear(void)
{
  for ( int ispin = 0; ispin < nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < nkp(); ikp++ )
    {
      if ( sd(ispin,ikp) != 0 )
      {
        if ( sdcontext(ispin,ikp)->active() )
        {
          sd(ispin,ikp)->c().clear();
        }
      }
    }
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
int Wavefunction::nempty() const { return nempty_; } // number of empty states

////////////////////////////////////////////////////////////////////////////////
int Wavefunction::nspin() const { return nspin_; } // number of empty states

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::entropy(void) const 
{
  assert(nspin_==1);
  assert(kpoint_.size()==1);
  return sd(0,0)->entropy(nspin_);
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
        if ( sdcontext_[ispin][ikp]!= 0 && sdcontext_[ispin][ikp]->active() )
        {
          sd_[ispin][ikp]->resize(cell,refcell,ecut,nst_[ispin]);
        }
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
void Wavefunction::reset(void)
{
  // reset all SlaterDets
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( sdcontext_[ispin][ikp]!= 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->reset();
      }
    }
  }
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
  reset();
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
  resize(cell_,refcell_,ecut_);
  reset();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nspin(int nspin)
{
  assert(nspin==1 || nspin==2);
  if ( nspin == nspin_ ) return;
  
  deallocate();
  cout << " Wavefunction::set_nspin: " << nspin << " deallocate done" << endl;
  
  nspin_ = nspin;
  
  compute_nst();
  allocate();
  cout << " Wavefunction::set_nspin: " << nspin << " allocate done" << endl;
  resize(cell_,refcell_,ecut_);
  reset();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nrowmax(int n)
{
  if ( n > ctxt_.size() )
  {
    cout << " Wavefunction::set_nrowmax: nrowmax > ctxt_.size()" << endl;
    return;
  }
  
  deallocate();
  nrowmax_ = n;
  compute_nst();
  allocate();
  resize(cell_,refcell_,ecut_);
  reset();
}
  
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::add_kpoint(D3vector kpoint, double weight)
{
  for ( int i = 0; i < kpoint_.size(); i++ )
  {
    if ( kpoint == kpoint_[i] )
    {
      cout << " Wavefunction::add_kpoint: warning: kpoint already defined" 
           << endl;
      //!! return;
    }
  }
  
  deallocate();
  cout << " Wavefunction::add_kpoint: " << kpoint << " deallocate done" << endl;
  
  kpoint_.push_back(kpoint);
  weight_.push_back(weight);
  
  allocate();
  cout << " Wavefunction::add_kpoint: " << kpoint << " allocate done" << endl;
  resize(cell_,refcell_,ecut_);
  reset();
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::del_kpoint(D3vector kpoint)
{
  cout << " Wavefunction::del_kpoint: not implemented" << endl;
  assert(false);
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::randomize(double amplitude)
{
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->randomize(amplitude);
      }
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
        if ( sd_[0][ikp] != 0 && sdcontext_[0][ikp]->active() )
        {
          sd_[0][ikp]->update_occ(nel_,nspin_);
        }
      }
    }
    else if ( nspin_ == 2 )
    {
      const int nocc_up = (nel_+1)/2+deltaspin_;
      const int nocc_dn = nel_/2 - deltaspin_;
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        if ( sd_[0][ikp] != 0 && sdcontext_[0][ikp]->active() )
          sd_[0][ikp]->update_occ(nocc_up,nspin_);
        if ( sd_[1][ikp] != 0 && sdcontext_[1][ikp]->active() )
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
    //!! double totalcharge = (double) ( nel_ + netcharge_ );
    const double totalcharge = (double) nel_;
    enum direction { up, down };
    direction dir = up;

    double rhosum = 0.0;
    for ( int ispin = 0; ispin < nspin_; ispin++ )
    {
      for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      {
        if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
        {
          sd_[ispin][ikp]->update_occ(nspin_,mu,temp);
          rhosum += sd_[ispin][ikp]->total_charge();
        }
        //!! rhosum must be reduced on pe 0 of each sdcontext only
        //!! without reduction, works only if nspin_==1 and nkp_==1
        assert(nspin_==1);
        assert(kpoint_.size()==1);
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
          if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
          {
            sd_[ispin][ikp]->update_occ(nspin_,mu,temp);
            rhosum += sd_[ispin][ikp]->total_charge();
          }
          //!! rhosum must be reduced on pe 0 of each sdcontext only
          //!! without reduction, works only if nspin_==1 and nkp_==1
          assert(nspin_==1);
          assert(kpoint_.size()==1);
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
      //!! print on one process only
      cout << " <!-- Wavefunction::update_occ: sum = "
           << rhosum << " -->" << endl;
      cout << " <!-- Wavefunction::update_occ: mu = "
           << setprecision(4) << mu / eVolt << " eV" << " -->" << endl;

      cout.setf(ios::right,ios::adjustfield);
      cout.setf(ios::fixed,ios::floatfield);
 
      cout << " <!-- Wavefunction::update_occ: occupation numbers" << endl;
      for ( int ispin = 0; ispin < nspin_; ispin++ )
      {
        for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
        {
          if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
          {
            for ( int n = 0; n < sd_[ispin][ikp]->nst(); n++ )
            {
              cout << setw(7) << setprecision(4) << sd_[ispin][ikp]->occ(n);
              if ( ( n%10 ) == 9 ) cout << endl;
            }
          }
        }
        cout << "  -->" << endl;
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
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->gram();
      }
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
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->riccati(*wf.sd_[ispin][ikp]);
      }
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
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->align(*wf.sd_[ispin][ikp]);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
double Wavefunction::dot(const Wavefunction& wf) const
{
  assert(wf.context() == ctxt_);
  double sum = 0.0;
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        sum += sd_[ispin][ikp]->dot(*wf.sd_[ispin][ikp]);
      }
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
    for ( int ikp = 0; ikp < nkp(); ikp++ )
    {
      if ( sd(ispin,ikp) != 0 )
      {
        if ( sdcontext(ispin,ikp)->active() )
        {
          // compute eigenvalues
          if ( sd(ispin,ikp)->basis().real() )
          {
            // proxy real matrices c, cp
            DoubleMatrix c(sd(ispin,ikp)->c());
            DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
            DoubleMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
            valarray<double> w(h.m());
 
            // factor 2.0 in next line: G and -G
            h.gemm('t','n',2.0,c,cp,0.0);
            // rank-1 update correction
            h.ger(-1.0,c,0,cp,0);
 
            // cout << " Hamiltonian at k = " << sd(ispin,ikp)->kpoint()
            //      << endl;
            // cout << h;
 
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
            if ( ctxt_.onpe0() )
            {
              const double eVolt = 2.0 * 13.6058;
              cout <<    "  <eigenvalues spin=\"" << ispin
                   << "\" kpoint=\"" << sd(ispin,ikp)->kpoint()
                   << "\" n=\"" << h.m() << "\">" << endl;
              for ( int i = 0; i < h.m(); i++ )
              {
                cout << setw(12) << setprecision(5) << w[i]*eVolt;
                if ( i%5 == 4 ) cout << endl;
              }
              if ( h.m()%5 != 0 ) cout << endl;
              cout << "  </eigenvalues>" << endl;
            }
            // set eigenvalues in SlaterDet
            sd(ispin,ikp)->set_eig(w);
          }
          else
          {
            // complex case not implemented
            assert(false);
            #if 0
            ComplexMatrix& c(wf.sd[ikp]->c());
            ComplexMatrix& cp(dwf.sd[ikp]->c());
 
            ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
            h.gemm('c','n',1.0,c,cp,0.0);
 
            //cout << " Hamiltonian at k = " << sd[ikp]->kpoint() << endl;
            //cout << h;
 
            valarray<double> w(h.m());

            h.heev('l',w);
            cout << " Eigenvalues at k = " << sd[ikp]->kpoint() << endl;
            const double eVolt = 2.0 * 13.6058;
            for ( int i = 0; i < h.m(); i++ )
            {
              cout << "%" << setw(3) << ikp
                   << setw(10) << setprecision(5) << w[i]*eVolt << endl;;
            }
            #endif
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::print(ostream& os, string encoding, string tag)
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
      sd_[ispin][ikp]->print(os,encoding);
  }
  
  if ( ctxt_.onpe0() )
    os << "</" << tag << ">" << endl;
}

#if USE_CSTDIO_LFS
////////////////////////////////////////////////////////////////////////////////
void Wavefunction::write(FILE* outfile, string encoding, string tag)
{
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
    off_t len = str.length();
    fwrite(str.c_str(),sizeof(char),len,outfile);
  }

  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      sd_[ispin][ikp]->write(outfile,encoding);
  }
  
  if ( ctxt_.onpe0() )
  {
    ostringstream os;
    os << "</" << tag << ">" << endl;
    string str(os.str());
    off_t len = str.length();
    fwrite(str.c_str(),sizeof(char),len,outfile);
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::info(ostream& os, string tag)
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
      sd_[ispin][ikp]->info(os);
  }
  
  if ( ctxt_.onpe0() )
    os << "</" << tag << ">" << endl;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, Wavefunction& wf)
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
  assert(cell_ == wf.cell_);
  assert(refcell_ == wf.refcell_);
  assert(ecut_ == wf.ecut_);
  
  weight_ = wf.weight_;
  kpoint_ = wf.kpoint_;
  
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( sd_[ispin][ikp] != 0 && sdcontext_[ispin][ikp]->active() )
      {
        *sd_[ispin][ikp] = *wf.sd_[ispin][ikp];
      }
    }
  }
  return *this;
}
