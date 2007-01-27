////////////////////////////////////////////////////////////////////////////////
//
// ConfinementPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConfinementPotential.h,v 1.2 2007-01-27 23:47:56 fgygi Exp $

#ifndef CONFINEMENTPOTENTIAL_H
#define CONFINEMENTPOTENTIAL_H

class Basis;
#include <valarray>
using namespace std;

class ConfinementPotential
{
  private:
  
  double ecuts_, facs_, sigmas_;
  const Basis& basis_;
  valarray<double> fstress_, dfstress_;

  public:
  
  double facs(void) const { return facs_; }
  double sigmas(void) const { return sigmas_; }
  double ecuts(void) const { return ecuts_; }
  
  const valarray<double>& fstress(void) const { return fstress_; }
  const valarray<double>& dfstress(void) const { return dfstress_; }
  
  void update(void);
  
  const Basis& basis() const { return basis_; }

  ConfinementPotential(double ecuts, double facs, double sigmas, 
    const Basis& basis);
  ~ConfinementPotential();
};
#endif
