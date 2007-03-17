////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.h,v 1.2 2007-03-17 01:14:00 fgygi Exp $

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Sample;
class EnergyFunctional;

#include <vector>
#include <valarray>

class Preconditioner
{
  private:
  
  const Sample& s_;
  const EnergyFunctional& ef_;
  std::vector<std::vector<std::valarray<double> > > diag_; // diag_[ispin][ikp][ig]

  public:

  void update(void);
  
  const std::valarray<double>& diag(int ispin, int ikp) const
  { return diag_[ispin][ikp]; }

  Preconditioner(const Sample& s, const EnergyFunctional& ef);
  //~Preconditioner();
};
#endif
