////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.C,v 1.5 2003-10-02 17:37:05 fgygi Exp $

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
  delete vcontext_;
  delete vbasis_;
  delete vft_;
  for ( int ikp = 0; ikp < ft_.size(); ikp++ )
    delete ft_[ikp];
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
      for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
      {
        if ( wf_.sd(ispin,ikp) != 0 && wf_.sdcontext(ispin,ikp)->active() )
        {
      
          wf_.sd(ispin,ikp)->compute_density(*ft_[ikp],
            wf_.weight(ikp), &rhor[ispin][0]);
 
        }
      }
    
      // sum across rows of spincontext[ispin] context
      //tmap["dsum"].start();
      wf_.spincontext(ispin)->dsum('r',vft_->np012loc(),1,
        &rhor[ispin][0],vft_->np012loc());
      //tmap["dsum"].stop();
      
      // check integral of charge density
      double sum = 0.0;
      for ( int i = 0; i < rhor[ispin].size(); i++ )
        sum += rhor[ispin][i];
      sum *= omega / vft_->np012();
      // cout << ctxt_.mype() << ": local integral of rhor: " << sum << endl;

      wf_.spincontext(ispin)->dsum('c',1,1,&sum,1);
      if ( ctxt_.onpe0() )
      {
        cout.setf(ios::fixed,ios::floatfield);
        cout.setf(ios::right,ios::adjustfield);
        cout << "  <total_electronic_charge> " << setprecision(8) << sum 
             << " </total_electronic_charge>" << endl;
      }

      // compute Fourier coefficients of the charge density
      for ( int i = 0; i < rhor[ispin].size(); i++ )
        rhotmp[i] = rhor[ispin][i];
      vft_->forward(&rhotmp[0],&rhog[ispin][0]);

      for ( int ig = 0; ig < vbasis_->localsize(); ig++ )
        rhog[ispin][ig] *= omega;

    }
  }
}

 
