////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.h,v 1.6 2007-03-17 01:14:00 fgygi Exp $

#ifndef NONLOCALPOTENTIAL_H
#define NONLOCALPOTENTIAL_H

#include "AtomSet.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "Context.h"
#include "Matrix.h"

class NonLocalPotential
{
  private:
  
  const Context& ctxt_;
  const AtomSet& atoms_;
  const SlaterDet& sd_;
  const Basis& basis_;

  int nsp;   // number of species
  int nspnl; // number of non-local species
  
  std::vector<int>                  lmax;     // lmax[is]
  std::vector<int>                  lloc;     // lloc[is]
  std::vector<int>                  na;       // na[is]
  std::vector<int>                  npr;      // npr[is]
  std::vector<int>                  nprna;    // nprna[is]
  std::vector<std::vector<int> >    lproj;    // lproj[is][ipr]
  std::vector<std::vector<double> > wt;       // wt[is][ipr]
  std::vector<std::vector<double> > twnl;     // twnl[is][npr*ngwl]
  std::vector<std::vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl],ij=0,..,5
  
  std::vector<int>             nquad;    // nquad[is]
  // rquad[is][iquad], iquad = 0, nquad[is]-1
  std::vector<std::vector<double> > rquad; 
  // wquad[is][iquad], iquad = 0, nquad[is]-1
  std::vector<std::vector<double> > wquad; 
    
  mutable TimerMap tmap;
  void init(void);
   
  public:
  
  NonLocalPotential(const AtomSet& as, const SlaterDet& sd) :  
    ctxt_(sd.context()), atoms_(as), sd_(sd), basis_(sd.basis()) { init(); }
  ~NonLocalPotential(void);
               
  void update_twnl(void);
  double energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, std::vector<std::vector<double> >& fion, 
    bool compute_stress, std::valarray<double>& sigma_enl);
};
#endif
