////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.h,v 1.2 2004-09-14 22:24:11 fgygi Exp $

#ifndef CHARGEDENSITY_H
#define CHARGEDENSITY_H

#include <vector>
#include <valarray>
#include <complex>
#include "Context.h"

class Wavefunction;
class FourierTransform;
class Basis;

class ChargeDensity
{
  private:
  
  const Context& ctxt_;
  const Wavefunction& wf_;
  const Context* vcontext_;
  Basis* vbasis_;
  FourierTransform* vft_;
  vector<FourierTransform*> ft_; // ft_[ikp];
  valarray<complex<double> > rhotmp;
  
  public:
  
  vector<vector<double> > rhor; // rhor[ispin][i]
  vector<vector<complex<double> > > rhog; // rhog[ispin][ig]

  void update_density(void);
  void update_rhor(void);
  
  Basis* vbasis(void) const { return vbasis_; }
  const Context* vcontext(void) const { return vcontext_; }
  FourierTransform* vft(void) const { return vft_; }
  FourierTransform* ft(int ikp) const { return ft_[ikp]; }

  ChargeDensity(const Wavefunction& wf);
  ~ChargeDensity();
};
#endif
