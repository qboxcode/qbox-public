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
// AndersonMixer.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <vector>
#include <valarray>
#include <cassert>
#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

class AndersonMixer
{
  // nmax is the dimension of the subspace of previous search directions
  // nmax=0: use simple mixing (no acceleration)
  // nmax=1: use one previous direction
  int     m_;                    // dimension of vectors
  int     nmax_;                 // maximum number of vectors (without current)
  int     n_;                    // number of vectors
  int     k_;                    // index of current vector
  const   MPI_Comm* const pcomm_;// pointer to relevant Context, null if local
  int     mype_;
  int     npes_;
  bool    diag_; // use diagonalization (default true)
  double  eig_ratio_; // eigenvalue ratio for regularization

  std::vector<std::valarray<double> > x_,f_;

  public:

  AndersonMixer(const int m, const int nmax, const MPI_Comm* const pcomm);
  void update(double* x, double* f, double* xbar, double* fbar);
  void restart(void);
  void set_diag(bool b) { diag_ = b; }
  void set_eig_ratio(double x) { eig_ratio_ = x; }
};
#endif
