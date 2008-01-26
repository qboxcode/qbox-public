////////////////////////////////////////////////////////////////////////////////
//
// SlaterDet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SlaterDet.h,v 1.22 2008-01-26 01:34:11 fgygi Exp $

#ifndef SLATERDET_H
#define SLATERDET_H

class FourierTransform;
#include "Context.h"
#include "Basis.h"
#include "Matrix.h"

#include "D3vector.h"
#include <vector>
#include <iostream>
#include "Timer.h"
#include <string>
#include <map>
#include "mpi.h"

typedef std::map<std::string,Timer> TimerMap;

class SlaterDet
{
  private:

  const Context& ctxt_;
  Context* my_col_ctxt_;
  Basis* basis_;
  ComplexMatrix c_;
  std::vector<double> occ_;
  std::vector<double> eig_;

  void byteswap_double(size_t n, double* x);
  double fermi(double e, double mu, double fermitemp);

  public:

  mutable TimerMap tmap;

  SlaterDet(const Context& ctxt, D3vector kpoint);
  SlaterDet(const SlaterDet& rhs);
  ~SlaterDet();
  const Context& context(void) const { return ctxt_; }
  const Basis& basis(void) const { return *basis_; }
  const D3vector kpoint(void) const { return basis_->kpoint(); }
  const ComplexMatrix& c(void) const { return c_; }
  ComplexMatrix& c(void) { return c_; }
  const std::vector<double>& occ(void) const { return occ_; }
  const std::vector<double>& eig(void) const { return eig_; }
  int nst(void) const { return c_.n(); }
  int nstloc(void) const { return c_.nloc(); }
  void resize(const UnitCell& cell, const UnitCell& refcell,
              double ecut, int nst);
  void compute_density(FourierTransform& ft, double weight, double* rho) const;
  void rs_mul_add(FourierTransform& ft, const double* v, SlaterDet& sdp) const;
  void randomize(double amplitude);
  void cleanup(void);
  void reset(void);
  void gram(void);
  void riccati(const SlaterDet& sd);
  void lowdin(void);
  void align(const SlaterDet& sd);
  void ortho_align(const SlaterDet& sd);
  std::complex<double> dot(const SlaterDet& sd) const;
  double total_charge(void);
  void update_occ(int nel, int nspin);
  void update_occ(int nspin, double mu, double temp);
  double eig(int i) const { return eig_[i]; };
  const double* eig_ptr(void) const { return &eig_[0]; }
  const double* eig_ptr(int i) const { return &eig_[i]; }
  double occ(int i) const { return occ_[i]; };
  const double* occ_ptr(void) const { return &occ_[0]; }
  const double* occ_ptr(int i) const { return &occ_[i]; }
  void set_occ(std::vector<double>& occ)
    { assert(occ_.size()==occ.size()); occ_ = occ; }
  void set_eig(std::vector<double>& eig)
    { assert(eig_.size()==eig.size()); eig_ = eig; }
  void set_eig(std::valarray<double>& eig)
    { assert(eig_.size()==eig.size());
      for ( int i = 0; i < eig.size(); i++ )
        eig_[i] = eig[i];
    }
  double entropy(int nspin);
  double ortho_error(void);
  double memsize(void) const;
  double localmemsize(void) const;
  SlaterDet& operator=(SlaterDet& rhs);
  void print(std::ostream& os, std::string encoding, double weight, int ispin,
    int nspin);
  void write(MPI_File& fh, std::string encoding, double weight, int ispin,
    int nspin);
  void info(std::ostream& os);
};
std::ostream& operator << ( std::ostream& os, SlaterDet& sd );

class SlaterDetException
{
  public:
  std::string msg;
  SlaterDetException(std::string s) : msg(s) {}
};
#endif
