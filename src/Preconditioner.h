////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.h,v 1.1 2004-03-11 21:58:10 fgygi Exp $

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

class Sample;
class EnergyFunctional;

#include <vector>
#include <valarray>
using namespace std;

class Preconditioner
{
  private:
  
  const Sample& s_;
  const EnergyFunctional& ef_;
  vector<vector<valarray<double> > > diag_; // diag_[ispin][ikp][ig]

  public:

  void update(void);
  
  const valarray<double>& diag(int ispin, int ikp) const
  { return diag_[ispin][ikp]; }

  Preconditioner(const Sample& s, const EnergyFunctional& ef);
  //~Preconditioner();
};
#endif
