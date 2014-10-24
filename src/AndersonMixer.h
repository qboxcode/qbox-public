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
// $Id: AndersonMixer.h,v 1.8 2009-03-08 01:10:30 fgygi Exp $

#ifndef ANDERSONMIXER_H
#define ANDERSONMIXER_H

#include <vector>
#include <valarray>
#include <cassert>
#include "Context.h"

class AndersonMixer
{
  // nmax is the dimension of the subspace of previous search directions
  // nmax=0: use simple mixing (no acceleration)
  // nmax=1: use one previous direction
  int     m_;                    // dimension of vectors
  int     nmax_;                 // maximum number of vectors (without current)
  int     n_;                    // number of vectors
  int     k_;                    // index of current vector
  const   Context* const pctxt_; // pointer to relevant Context, null if local

  std::vector<std::valarray<double> > x_,f_;

  public:

  AndersonMixer(const int m, const int nmax, const Context* const pctxt) :
    m_(m), nmax_(nmax), pctxt_(pctxt)
  {
    assert( nmax >= 0 );
    x_.resize(nmax_+1);
    f_.resize(nmax_+1);
    for ( int n = 0; n < nmax_+1; n++ )
    {
      x_[n].resize(m_);
      f_[n].resize(m_);
    }
    restart();
  }

  void update(double* x, double* f, double* xbar, double* fbar);
  void restart(void);
};
#endif
