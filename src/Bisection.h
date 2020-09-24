////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011 The Regents of the University of California
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
// Bisection.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BISECTION_H
#define BISECTION_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>

#include "Context.h"
#include "SlaterDet.h"
#include "Matrix.h"

class FourierTransform;
class Bisection
{
  private:

    const Context& ctxt_;

    // bisection levels in each directions
    int nlevels_[3]; // bisection level
    int ndiv_[3];    // number of subdivisions
    int nlevelsmax_; // max level of bisection
    int nst_;        // number of states in SlaterDet
    int nstloc_;

    // real space grid
    int np_[3];            // grid size
    int np01_;             // size of one xy plane
    int np2loc_;           // local z size
    int np012loc_;         // local size
    FourierTransform *ft_;

    std::vector<long int> localization_; // localization indices

    // xy_proj[ir]: projector in xy plane associated with grid point ir
    std::vector<int> xy_proj_;

    // matrices of real space wave functions in subdomains
    std::vector<DoubleMatrix*> rmat_;

    // a matrices
    int nmat_;
    std::vector<DoubleMatrix*> amat_;
    std::vector<std::vector<double> > adiag_;
    DoubleMatrix *u_;

    // test function
    bool check_amat(const ComplexMatrix &c);
    void trim_amat(const std::vector<double>& occ);

  public:

    Bisection(const SlaterDet& sd, const int nlevels[3]);
    void compute_transform(const SlaterDet& sd, int maxsweep, double tol);
    void compute_localization(double epsilon);
    void forward(SlaterDet& sd);
    void forward(DoubleMatrix& u, SlaterDet& sd);
    void backward(SlaterDet& sd);
    void backward(DoubleMatrix& u, SlaterDet& sd);

    int nmat(void) const { return nmat_; }
    long int localization(int i) const { return localization_[i]; }
    const std::vector<long int>& localization(void) const
    { return localization_; }
    bool overlap(int i, int j) const;
    bool overlap(const std::vector<long int>& loc, int i, int j) const;
    const DoubleMatrix& u(void) const { return *u_; }
    double pair_fraction(void) const;
    double size(int i) const;
    double total_size(void) const;
    ~Bisection();
};
#endif
