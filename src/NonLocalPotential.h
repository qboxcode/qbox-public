////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.h,v 1.3 2004-02-04 19:55:17 fgygi Exp $

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
  
  vector<int>             lmax;     // lmax[is]
  vector<int>             lloc;     // lloc[is]
  vector<int>             na;       // na[is]
  vector<int>             naloc;    // naloc[is]
  vector<int>             nalocmax; // nalocmax[is]
  vector<int>             npr;      // npr[is]
  vector<int>             nprna;    // nprna[is]
  vector<vector<int> >    lproj;    // lproj[is][ipr]
  vector<vector<double> > wt;       // wt[is][ipr]
  vector<vector<double> > twnl;     // twnl[is][npr*ngwl]
  vector<vector<double> > dtwnl;    // dtwnl[is][6*npr*ngwl], ij=0,..,5
  
  vector<DoubleMatrix*>   anl;      // anl[is][ipr*ia][ig]
  vector<vector<double> > singr;    // singr[is][naloc*2*ngwloc]
  vector<vector<double> > cosgr;    // cosgr[is][naloc*2*ngwloc]
  
  vector<int>             nquad;    // nquad[is]
  vector<vector<double> > rquad;    // rquad[is][iquad], iquad = 0, nquad[is]-1
  vector<vector<double> > wquad;    // wquad[is][iquad], iquad = 0, nquad[is]-1
    
  void init(void);
   
  public:
  
  NonLocalPotential(const AtomSet& as, const SlaterDet& sd) :  
    ctxt_(sd.context()), atoms_(as), sd_(sd), basis_(sd.basis()) { init(); }
  ~NonLocalPotential(void);
               
  void update_twnl(void);
  void update_eigr(vector<vector<double> >& tau);
  void update_anl(void);
  double energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl);
};
#endif
