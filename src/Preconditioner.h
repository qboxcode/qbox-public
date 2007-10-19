////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.h,v 1.3 2007-10-19 16:24:04 fgygi Exp $

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
