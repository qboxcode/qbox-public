////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.C,v 1.10 2004-10-28 16:48:15 fgygi Exp $

#include "ChargeDensity.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "SlaterDet.h"

#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::ChargeDensity(const Wavefunction& wf) : ctxt_(wf.context()),
wf_(wf) 
{
  // determine which k points are stored on the sd context to which 
  // this pe belongs
  int myspin = -1;
  int ikplast = -1;
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      if ( wf.sd(ispin,ikp) != 0 )
      {
        myspin = ispin;
        ikplast = ikp;
      }
    }
  }
  // the last kpoint on this process is ikplast
  
  // create a vbasis
  // Note: avoid creating more than one vcontext if there are 
  // more than one kpoint per process
  if ( wf.sd(myspin,ikplast) != 0 )
  {
    // Note: Basis defined on columns of vcontext
    vcontext_ = new Context(*wf_.spincontext(myspin),'c',
      wf_.spincontext(myspin)->mycol());
    vbasis_ = new Basis(*vcontext_, D3vector(0,0,0));
    vbasis_->resize(wf.cell(),wf.refcell(),4.0*wf.ecut());
    //cout << vcontext_->mype() << ": vcontext = " << *vcontext_ << endl;
  }
  else
  {
    vcontext_ = 0;
    vbasis_ = 0;
  }
    
  // define vft_, FT on vbasis context for transforming the density
  
  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
  vft_ = new FourierTransform(*vbasis_,np0v,np1v,np2v);
  
  rhor.resize(wf.nspin());
  rhog.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    rhor[ispin].resize(vft_->np012loc());
    rhog[ispin].resize(vbasis_->localsize());
  }
  rhotmp.resize(vft_->np012loc());
  
  // FT for interpolation of wavefunctions on the fine grid
  ft_.resize(wf.nkp());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      if ( wf.sd(ispin,ikp) != 0 )
        ft_[ikp] =
          new FourierTransform(wf.sd(ispin,ikp)->basis(),np0v,np1v,np2v);
      else
        ft_[ikp] = 0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::~ChargeDensity(void)
{
  delete const_cast<Context*>(vcontext_);
  delete vbasis_;
  delete vft_;
  for ( int ikp = 0; ikp < ft_.size(); ikp++ )
    delete ft_[ikp];
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    if ( ctxt_.myproc()==0 )
    {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_density(void)
{
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();
  
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if ( wf_.spincontext(ispin)->active())
    {
      assert(rhor[ispin].size() == vft_->np012loc() );
      assert(rhotmp.size() == vft_->np012loc() );
      tmap["charge_compute"].start();
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if ( wf_.sd(ispin,ikp) != 0 && wf_.sdcontext(ispin,ikp)->active() )
        {
      
          wf_.sd(ispin,ikp)->compute_density(*ft_[ikp],
            wf_.weight(ikp), &rhor[ispin][0]);
 
        }
      }
      tmap["charge_compute"].stop();
    
      // sum across rows of spincontext[ispin] context
      tmap["charge_rowsum"].start();
      wf_.spincontext(ispin)->dsum('r',vft_->np012loc(),1,
        &rhor[ispin][0],vft_->np012loc());
      tmap["charge_rowsum"].stop();
      
      // check integral of charge density
      // compute Fourier coefficients of the charge density
      double sum = 0.0;
      const int rhor_size = rhor[ispin].size();
      const double *const prhor = &rhor[ispin][0];
      tmap["charge_integral"].start();
      #pragma ivdep
      for ( int i = 0; i < rhor_size; i++ )
      {
        const double prh = prhor[i];
        sum += prh;
        rhotmp[i] = complex<double>(omega * prh, 0.0);
      }
      sum *= omega / vft_->np012();

      wf_.spincontext(ispin)->dsum('c',1,1,&sum,1);
      tmap["charge_integral"].stop();
      if ( ctxt_.onpe0() )
      {
        cout.setf(ios::fixed,ios::floatfield);
        cout.setf(ios::right,ios::adjustfield);
        cout << "  <!-- total_electronic_charge: " << setprecision(8) << sum 
             << " -->" << endl;
      }
        
      tmap["charge_ft"].start();
      vft_->forward(&rhotmp[0],&rhog[ispin][0]);
      tmap["charge_ft"].stop();
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_rhor(void)
{
  // recalculate rhor from rhog
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();
  assert(omega!=0.0);
  const double omega_inv = 1.0 / omega;
  
  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    if ( wf_.spincontext(ispin)->active())
    {
      assert(rhor[ispin].size() == vft_->np012loc() );
      assert(rhotmp.size() == vft_->np012loc() );
      
      vft_->backward(&rhog[ispin][0],&rhotmp[0]);

      const int rhor_size = rhor[ispin].size();
      double *const prhor = &rhor[ispin][0];
      #pragma ivdep
      for ( int i = 0; i < rhor_size; i++ )
      {
        prhor[i] = rhotmp[i].real() * omega_inv;
      }
    }
  }
}

 
