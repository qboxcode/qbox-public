////////////////////////////////////////////////////////////////////////////////
//
// Wavefunction.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Wavefunction.C,v 1.11 2003-10-02 17:33:20 fgygi Exp $

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
deltaspin_(wf.deltaspin_), cell_(wf.cell_), refcell_(wf.refcell_),
ecut_(wf.ecut_), weight_(wf.weight_), kpoint_(wf.kpoint_)
{
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
  
  resize();
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
  
  int npr, npc; // dimensions of each sdcontext
  npr = ctxt_.size() / nkp;
  npc = 1;
  while ( npr%2 == 0 && npr > nrowmax_ )
  {
    npr /= 2;
    npc *= 2;
  }
  spincontext_.resize(1);
  sdcontext_.resize(1);
  sdcontext_[0].resize(1);
  sd_.resize(1);
  sd_[0].resize(1);
  
  spincontext_[0] = new Context(ctxt_,npr,npc);
  sdcontext_[0][0] = 0;
  sd_[0][0] = 0;
  if ( spincontext_[0]->active() )
  {
    sdcontext_[0][0] = new Context(*spincontext_[0]);
    sd_[0][0] = new SlaterDet(*sdcontext_[0][0],kpoint_[0]);
  }
  
#if 0
  cout << ctxt_.mype() << ": Wavefunction::allocate: start " 
       << "nspin = " << nspin_ << " nkp = " << nkp << endl;
  ctxt_.barrier();
  
  spincontext_.resize(nspin_);
  
  int npr, npc; // dimensions of each sdcontext
  
  // create spincontext_[ispin]
  if ( nspin_ == 1 )
  {
    // determine dimensions of spincontext
    if ( ctxt_.size() >= nkp )
    {
      // more than one task per k point
      npr = ctxt_.size() / nkp;
      npc = 1;
      while ( npr%2 == 0 && npr > nrowmax_ )
      {
        npr /= 2;
        npc *= 2;
      }
      spincontext_[0] = new Context(ctxt_,npr,nkp*npc);
    }
    else
    {
      // less than one task per k point
      npr = 1;
      npc = 1;
      spincontext_[0] = new Context(ctxt_,'r');
    }
  }
  else
  {
    // nspin == 2
    nst_[0] = ( nel_ + 1 ) / 2 + deltaspin_ + nempty_;
    nst_[1] = nel_ / 2 - deltaspin_ + nempty_;
    
    if ( ctxt_.size() > 1 )
    {
      // determine dimensions of subctxt
      if ( ctxt_.size() >= 2*nkp )
      {
        // more than one task per k point
        npr = ctxt_.size() / (2*nkp);
        npc = 1;
//         if ( npr%2 == 0 && npr > nrowmax )
//         {
//           npr /= 2;
//           npc *= 2;
//         }
        spincontext_[0] = new Context(ctxt_,0,npr,nkp*npc);
        spincontext_[1] = new Context(ctxt_,nkp*npr*npc,npr,nkp*npc);
      }
      else
      {
        // less than one task per k point
        npr = 1;
        npc = 1;
        Context subctxt(ctxt_,2,ctxt_.size()/2);
        if ( subctxt.active() )
        {
          spincontext_[0] = new Context(subctxt,'r',0);
          spincontext_[1] = new Context(subctxt,'r',1);
        }
        else
        {
          spincontext_[0] = 0;
          spincontext_[1] = 0;
        }
      }
    }
    else
    {
      // only 1 task. both spins reside on the same node
      // both spincontexts are copies of ctxt_
      npr = 1;
      npc = 1;
      spincontext_[0] = new Context(ctxt_);
      spincontext_[1] = new Context(ctxt_);
    }
  }
  ctxt_.barrier();
  
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    if ( spincontext_[ispin] != 0 )
    {
      cout << ctxt_.mype() << ": spincontext[" << ispin << "]: "
           << *spincontext_[ispin];
    }
  }
  ctxt_.barrier();

  sdcontext_.resize(nspin_);
  sd_.resize(nspin_);
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    sdcontext_[ispin].resize(nkp);
    sd_[ispin].resize(nkp);
  }
    
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    if ( spincontext_[ispin] != 0 && spincontext_[ispin]->active() )
    {
      if ( spincontext_[ispin]->size() >= nkp )
      {
        for ( int ikp = 0; ikp < nkp; ikp++ )
        {
          ctxt_.barrier();
          cout << ctxt_.mype()
            << ": Wavefunction::allocate: creating sdcontext "
            << " ispin=" << ispin 
            << " ikp=" << ikp << " npr=" << npr << " npc=" << npc << endl;
//           sdcontext_[ispin][ikp] = 
//             new Context(*spincontext_[ispin],ikp*npr*npc,npr,npc);
// !! note: next line works only for nkp == 1
          sdcontext_[ispin][ikp] = 
            new Context(*spincontext_[ispin],0,ikp*npc,npr,npc,npr,npc);
          if ( sdcontext_[ispin][ikp]->active() )
            sd_[ispin][ikp] =
              new SlaterDet(*sdcontext_[ispin][ikp],kpoint_[ikp]);
          else
            sd_[ispin][ikp] = 0;
        }
      }
      else
      {
        // more than one k point per task
        // round robin allocation of ikp to contexts
        int icol = 0;
        const int npcol = spincontext_[ispin]->npcol(); 
        for ( int ikp = 0; ikp < nkp; ikp++ )
        {
          assert( icol < npcol );
          if ( spincontext_[ispin]->active() )
          {
            sdcontext_[ispin][ikp] = new Context(*spincontext_[ispin],'c',icol);
            if ( sdcontext_[ispin][ikp]->active() )
              sd_[ispin][ikp] =
                new SlaterDet(*sdcontext_[ispin][ikp],kpoint_[ikp]);
          }
          else
          {
            sdcontext_[ispin][ikp] = 0;
            sd_[ispin][ikp] = 0;
          }
          // loop back to 0 if past npc-1
          icol = ( icol + 1 ) % npcol;
        }
      }
      cout << ctxt_.mype() << ": Wavefunction::allocate: done"
           << " ispin = " << ispin << " nkp = " << nkp << endl;
    }
    else
    {
      // *spincontext_[ispin] is inactive
      for ( int ikp = 0; ikp < nkp; ikp++ )
      {
        sdcontext_[ispin][ikp] = 0;
        sd_[ispin][ikp] = 0;
      }
    }
  } // ispin
  cout << ctxt_.mype() << ": Wavefunction::allocate: end " 
       << "nspin = " << nspin_ << " nkp = " << nkp << endl;
#endif

#if 0
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < nkp; ikp++ )
    {
      if ( sdcontext_[ispin][ikp] != 0 )
      {
        cout << "<!-- ";
        cout << ctxt_.mype() << ": sdcontext[" << ispin << "][" << ikp << "]: "
           << *sdcontext_[ispin][ikp];
        cout << " -->" << endl;
      }
    }
  }
#endif
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
double Wavefunction::entropy(void) const {return 0.0;}// dimensionless entropy

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::resize(const UnitCell& cell, const UnitCell& refcell, 
  double ecut)
{
  cell_ = cell;
  refcell_ = refcell;
  ecut_ = ecut;
  
  resize();
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
void Wavefunction::resize(void)
{
  // resize SlaterDets if ecut,cell,refcell,or nst have changed
  
  // resize all SlaterDets using cell_, refcell_, ecut_ and nst_[ispin]
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
    {
      if ( sdcontext_[ispin][ikp]!= 0 && sdcontext_[ispin][ikp]->active() )
      {
        sd_[ispin][ikp]->resize(cell_,refcell_,ecut_,nst_[ispin]);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::set_nel(int nel)
{
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
  if ( nempty < 0 )
  {
    cout << " Wavefunction::set_nempty: negative value" << endl;
    return;
  }
  nempty_ = nempty;
  compute_nst();
  resize();
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
  resize();
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
  allocate();
  resize();
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
  resize();
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
  // update occupation numbers at finite temperature temp
  const double eVolt = 0.036749023; // 1 eV in Hartree
  
  // loop to find value of mu
  double mu = 0.0;
  double dmu = 1.0 * eVolt;
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
  while ( niter < 50 && fabs(rhosum - nel_) > 1.e-10 )
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
  
  if ( niter == 50 )
  {
    cout << "Wavefunction::update_occ: mu did not converge in 50 iterations"
         << endl;
    ctxt_.abort(1);
  }
  
  cout << " Wavefunction::update_occ: sum = " << rhosum << endl;

  cout << " Wavefunction::update_occ: mu = "
       << setprecision(4) << mu / eVolt << " eV" << endl;

  cout.setf(ios::right,ios::adjustfield);
  cout.setf(ios::fixed,ios::floatfield);
  
  //!! print on all processes
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
    cout << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Wavefunction::update_occ(void)
{
  // update occupation numbers for zero temperature
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
void Wavefunction::print(ostream& os, string encoding, string tag)
{
  os << "<" << tag << " ecut=\"" << ecut_ << "\""
     << " nspin=\"" << nspin_ << "\""
     << " nel=\"" << nel_ << "\""
     << " nempty=\"" << nempty_ << "\">" << endl;
  os << "<domain a=\"" 
     << cell_.a(0) << "\"\n        b=\""
     << cell_.a(1) << "\"\n        c=\""
     << cell_.a(2) << "\"/>" << endl;
  os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
     <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
     <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
     
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      sd_[ispin][ikp]->print(os,encoding);
  }
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
    os << "<domain a=\""
       << cell_.a(0) << "\"\n        b=\""
       << cell_.a(1) << "\"\n        c=\""
       << cell_.a(2) << "\"/>" << endl;
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
  os << "<" << tag << " ecut=\"" << ecut_ << "\""
     << " nspin=\"" << nspin_ << "\""
     << " nel=\"" << nel_ << "\""
     << " nempty=\"" << nempty_ << "\">" << endl;
  os << "<domain a=\"" 
     << cell_.a(0) << "\"\n        b=\""
     << cell_.a(1) << "\"\n        c=\""
     << cell_.a(2) << "\"/>" << endl;
  os << "<grid nx=\"" << sd_[0][0]->basis().np(0) << "\""
     <<      " ny=\"" << sd_[0][0]->basis().np(1) << "\""
     <<      " nz=\"" << sd_[0][0]->basis().np(2) << "\"/>" << endl;
     
  for ( int ispin = 0; ispin < nspin_; ispin++ )
  {
    for ( int ikp = 0; ikp < kpoint_.size(); ikp++ )
      sd_[ispin][ikp]->info(os);
  }
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
