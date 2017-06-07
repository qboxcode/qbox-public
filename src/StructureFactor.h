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
// StructureFactor.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STRUCTUREFACTOR_H
#define STRUCTUREFACTOR_H

#include <vector>
#include <complex>
#include <iostream>

class Basis;

class StructureFactor
{
  private:

  int _nsp, _ng;
  std::vector<int> _na;

  int _k0max, _k1max, _k2max,
      _k0min, _k1min, _k2min,
      _k0range, _k1range, _k2range;

  public:

  // convenience pointer access functions:
  // double *c0 = cos0_ptr(is,ia);
  // c0[ kx ] == cos(-i gx*tau[is][ia].x)

  double *cos0_ptr(int is, int ia)
    { return &cos0[is][ia*_k0range-_k0min]; }

  double *cos1_ptr(int is, int ia)
    { return &cos1[is][ia*_k1range-_k1min]; }

  double *cos2_ptr(int is, int ia)
    { return &cos2[is][ia*_k2range-_k2min]; }

  double *sin0_ptr(int is, int ia)
    { return &sin0[is][ia*_k0range-_k0min]; }

  double *sin1_ptr(int is, int ia)
    { return &sin1[is][ia*_k1range-_k1min]; }

  double *sin2_ptr(int is, int ia)
    { return &sin2[is][ia*_k2range-_k2min]; }

  // kx in [k0min, k0max]
  // ky in [k1min, k1max]
  // kz in [k2min, k2max]

  std::vector<std::vector<double> > cos0;  // cos0[is][ia*k0range-k0min+kx]
  std::vector<std::vector<double> > cos1;  // cos1[is][ia*k1range-k1min+ky]
  std::vector<std::vector<double> > cos2;  // cos2[is][ia*k2range-k2min+ky]
  std::vector<std::vector<double> > sin0;  // sin0[is][ia*k0range-k0min+kx]
  std::vector<std::vector<double> > sin1;  // sin1[is][ia*k1range-k1min+ky]
  std::vector<std::vector<double> > sin2;  // sin2[is][ia*k2range-k2min+ky]
  std::vector<std::vector<std::complex<double> > > sfac;  // sfac[is][ig]

  void init(const std::vector<std::vector<double> >& tau, const Basis& basis);
  void update(const std::vector<std::vector<double> >& tau, const Basis& basis);

};
#endif
