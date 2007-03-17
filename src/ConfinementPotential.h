////////////////////////////////////////////////////////////////////////////////
//
// ConfinementPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConfinementPotential.h,v 1.3 2007-03-17 01:14:00 fgygi Exp $

#ifndef CONFINEMENTPOTENTIAL_H
#define CONFINEMENTPOTENTIAL_H

#include <valarray>
class Basis;

class ConfinementPotential
{
  private:
  
  double ecuts_, facs_, sigmas_;
  const Basis& basis_;
  std::valarray<double> fstress_, dfstress_;

  public:
  
  double facs(void) const { return facs_; }
  double sigmas(void) const { return sigmas_; }
  double ecuts(void) const { return ecuts_; }
  
  const std::valarray<double>& fstress(void) const { return fstress_; }
  const std::valarray<double>& dfstress(void) const { return dfstress_; }
  
  void update(void);
  
  const Basis& basis() const { return basis_; }

  ConfinementPotential(double ecuts, double facs, double sigmas, 
    const Basis& basis);
  ~ConfinementPotential();
};
#endif
