////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SlaterDet.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SLATERDET_H
#define SLATERDET_H

class FourierTransform;
#include "Context.h"
#include "Basis.h"
#include "Matrix.h"

#include "D3vector.h"
#include <iosfwd>
#include "Timer.h"
#include <string>
#include <map>

class SharedFilePtr;

typedef std::map<std::string,Timer> TimerMap;

class SlaterDet
{
  private:

  const Context& ctxt_;
  MPI_Comm my_col_comm_;
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
  void compute_tau(FourierTransform& ft, double weight, double* taur) const;
  void rs_mul_add(FourierTransform& ft, const double* v, SlaterDet& sdp) const;
  void randomize(double amplitude);
  void cleanup(void);
  void init(void);
  void gram(void);
  void riccati(const SlaterDet& sd);
  void lowdin(void);
  void align(const SlaterDet& sd);
  void ortho_align(const SlaterDet& sd);
  std::complex<double> dot(const SlaterDet& sd) const;
  double total_charge(void) const;
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
  void set_occ(int i, double f) { occ_[i] = f; }
  void set_eig(std::vector<double>& eig)
    { assert(eig_.size()==eig.size()); eig_ = eig; }
  void set_eig(std::valarray<double>& eig)
    { assert(eig_.size()==eig.size());
      for ( int i = 0; i < eig.size(); i++ )
        eig_[i] = eig[i];
    }
  double entropy(int nspin) const;
  double ortho_error(void) const;
  double memsize(void) const;
  double localmemsize(void) const;
  SlaterDet& operator=(SlaterDet& rhs);
  void print(std::ostream& os, std::string encoding, double weight, int ispin,
    int nspin) const;
  void write(SharedFilePtr& fh, std::string encoding, double weight, int ispin,
    int nspin) const;
  void info(std::ostream& os) const;
  double empty_row_error(void);
  double g0_imag_error(void);
};
std::ostream& operator << ( std::ostream& os, SlaterDet& sd );

class SlaterDetException
{
  public:
  std::string msg;
  SlaterDetException(std::string s) : msg(s) {}
};
#endif
