////////////////////////////////////////////////////////////////////////////////
//
// ChargeDensity.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.h,v 1.6 2007-10-19 16:24:04 fgygi Exp $

#ifndef CHARGEDENSITY_H
#define CHARGEDENSITY_H

#include <vector>
#include <valarray>
#include <complex>
#include <string>
#include <map>
#include "Timer.h"
#include "Context.h"

class Wavefunction;
class FourierTransform;
class Basis;

typedef std::map<std::string,Timer> TimerMap;

class ChargeDensity
{
  private:

  const Context& ctxt_;
  const Context& vcontext_;
  const Wavefunction& wf_;
  Basis* vbasis_;
  FourierTransform* vft_;
  std::vector<FourierTransform*> ft_; // ft_[ikp];
  std::valarray<std::complex<double> > rhotmp;

  public:

  mutable TimerMap tmap;

  std::vector<std::vector<double> > rhor; // rhor[ispin][i]
  std::vector<std::vector<std::complex<double> > > rhog; // rhog[ispin][ig]

  void update_density(void);
  void update_rhor(void);

  Basis* vbasis(void) const { return vbasis_; }
  const Context& vcontext(void) const { return vcontext_; }
  FourierTransform* vft(void) const { return vft_; }
  FourierTransform* ft(int ikp) const { return ft_[ikp]; }

  ChargeDensity(const Wavefunction& wf);
  ~ChargeDensity();
};
#endif
