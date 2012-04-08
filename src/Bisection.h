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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>

#include "Context.h"
#include "Sample.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "SlaterDet.h"
#include "Matrix.h"
#include "Timer.h"
#include "isodate.h"
#include "jade.h"

using namespace std;

class Bisection
{
  private:

    Sample& s_;
    Context gcontext_;

    // bisection levels in each directions
    int nlevels_[3]; // bisection level
    int ndiv_[3];    // number of subdivisions
    int nlevelsmax_; // max level of bisection

    // real space grid
    int np_[3];            // grid size
    int np01_;             // size of one xy plane
    int np2loc_;           // local z size
    int np012loc_;         // local size
    FourierTransform *ft_;

    vector<vector<long int> > localization_;

    // xy_proj[ir]: projector in xy plane associated with grid point ir
    vector<int> xy_proj_;

    // matrices of real space wave functions in subdomains
    vector<DoubleMatrix*> rmat_;

    // a matrices
    int nmat_;
    vector<DoubleMatrix*> amat_;
    DoubleMatrix *u_;

    // test function
    bool check_amat(ComplexMatrix &c);
    void trim_amat(const vector<double>& occ);
    void distribute(int ispin);

  public:

    Bisection(Sample& s, int nlevels[3]);
    int nmat(void) const { return nmat_; }
    void localize(Wavefunction &wf, double epsilon);
    long int localization(int ispin, int i) const
    { return localization_[ispin][i]; }
    vector<vector<long int> > localization(void) const { return localization_; }
    bool overlap(int ispin, int i, int j);
    double pair_fraction(int ispin);
    double size(int ispin, int i);
    double total_size(int ispin);
    ~Bisection();
};
