////////////////////////////////////////////////////////////////////////////////
//
// Wavefunction.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Wavefunction.h,v 1.16 2007-10-19 16:24:05 fgygi Exp $

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "D3vector.h"
#include "UnitCell.h"
#include <vector>
#if USE_CSTDIO_LFS
#include <cstdio>
#endif

class SlaterDet;
class Context;

class Wavefunction
{
  private:

  const Context& ctxt_;

  int nel_;           // number of electrons
  int nempty_;        // number of empty states
  int nspin_;         // number of spins (1 or 2)
  int deltaspin_;     // number of spin excitations

  int nrowmax_;       // maximum number of rows of a spincontext

  UnitCell cell_ ;    // unit cell
  UnitCell refcell_ ; // reference cell
  double   ecut_ ;    // energy cutoff

  std::vector<double>    weight_;  // weight[ikp]
  std::vector<D3vector>  kpoint_;  // kpoint[ikp]

  std::vector<int> nst_;                       // nst_[ispin]
  std::vector<Context*> spincontext_;          // spincontext[ispin]
  std::vector<std::vector<Context*> > sdcontext_;   // sdcontext_[ispin][ikp]
  std::vector<std::vector<SlaterDet*> > sd_;        // sd[ispin][ikp]

  void allocate(); // create contexts and allocate SlaterDet's
  void deallocate();
  void compute_nst();
  void resize(); // resize SlaterDets if ecut,cell,refcell,or nst have changed

  public:

  Wavefunction(const Context& ctxt);
  Wavefunction(const Wavefunction& wf);
  ~Wavefunction();
  Wavefunction& operator=(const Wavefunction& wf);

  const Context& context(void) const { return ctxt_; }
  const UnitCell& cell(void) const { return cell_; }
  const UnitCell& refcell(void) const { return refcell_; }
  const D3vector kpoint(int ikp) const { return kpoint_[ikp]; }
  double weight(int ikp) const { return weight_[ikp]; }
  double ecut(void) const { return ecut_; }
  SlaterDet* sd(int ispin, int ikp) const { return sd_[ispin][ikp]; }
  const Context* sdcontext(int ispin, int ikp) const
    { return sdcontext_[ispin][ikp]; }
  const Context* spincontext(int ispin) const
    { return spincontext_[ispin]; }
  int nkp(void) const;            // number of k points
  int nel(void) const;            // total number of electrons
  int nst(int ispin) const;       // number of states of spin ispin
  int nst(void) const;            // number of states
  int nempty(void) const;         // number of empty states
  int nspin(void) const;          // number of spins
  int deltaspin(void) const;      // number of spin excitations
  int nrowmax(void) const { return nrowmax_; }

  double spin(void) const;        // total spin

  void resize(const UnitCell& cell, const UnitCell& refcell, double ecut);
  void resize(double ec) { resize(cell_,refcell_,ec); }
  void reset(void); // initialize with lowest plane waves
  void clear(void); // initialize with zero
  void set_nel(int nel);
  void set_nempty(int nempty);
  void set_nspin(int nspin);
  void set_deltaspin(int deltaspin);
  void set_nrowmax(int n);
  void add_kpoint(D3vector kpoint, double weight);
  void del_kpoint(D3vector kpoint);

  void randomize(double amplitude);

  void update_occ(double fermitemp);
  double entropy(void) const; // dimensionless entropy
  void gram(void);
  void riccati(Wavefunction& wf);
  void align(Wavefunction& wf);
  void diag(Wavefunction& dwf, bool eigvec);

  double dot(const Wavefunction& wf) const;

  void print(std::ostream& os, std::string encoding, std::string tag) const;
#if USE_CSTDIO_LFS
  void write(FILE* outfile, std::string encoding, std::string tag) const;
#endif
  void info(std::ostream& os, std::string tag) const;
};
std::ostream& operator << ( std::ostream& os, const Wavefunction& wf );
#endif
