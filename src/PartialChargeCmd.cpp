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
// PartialChargeCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PartialChargeCmd.h"

#include <iostream>
#include <complex>

#include "Context.h"
#include "Sample.h"
#include "Basis.h"
#include "Atom.h"
#include "ChargeDensity.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int PartialChargeCmd::action(int argc, char **argv)
{
  string usage("  Use: partial_charge [-spin {1|2}] name radius");

  // parse arguments
  // plot [-spin {1|2}] name radius

  // ispin = 0: include both spins
  // ispin = 1: include first spin
  // ispin = 2: include second spin
  int ispin = 0;
  double radius = 0.0;
  string atom_name;

  if ( !(argc==3 || argc==5) )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }

  const int nspin = s->wf.nspin();
  // process arguments
  int iarg = 1;
  if ( !strcmp(argv[iarg],"-spin") )
  {
    if ( nspin != 2 )
    {
      if ( ui->onpe0() )
        cout << "nspin = 1, cannot select spin" << endl;
      return 1;
    }
    iarg++;
    // process argument: ispin
    if ( iarg==argc )
    {
      if ( ui->onpe0() )
        cout << usage << endl;
      return 1;
    }
    ispin = atoi(argv[iarg]);
    if ( ispin < 1 || ispin > 2 )
    {
      if ( ui->onpe0() )
        cout << " spin must be 1 or 2" <<  endl;
      return 1;
    }
    iarg++;
  }

  // argument must be the atom name followed by the radius
  if ( iarg==argc )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }
  atom_name = argv[iarg];
  iarg++;
  if ( iarg==argc )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }
  radius = atof(argv[iarg]);

  // find atom and check validity of radius argument
  Atom* a = s->atoms.findAtom(atom_name);
  if ( a == 0 )
  {
    if ( ui->onpe0() )
      cout << " PartialChargeCmd: atom " << atom_name << " not defined" << endl;
    return 1;
  }
  if ( radius <= 0.0 )
  {
    if ( ui->onpe0() )
      cout << " PartialChargeCmd: radius must be positive" << endl;
    return 1;
  }

  D3vector pos = a->position();
  if ( ui->onpe0() )
  {
    cout << setprecision(5) << " Atom " << atom_name << " at " << pos << endl;
    cout << " radius = " << radius << endl;
  }

  const Context& ctxt = s->wf.sd_context();
  ChargeDensity cd(s->wf);
  Basis *vbasis = cd.vbasis();
  cd.update_density();
  const double omega = vbasis->cell().volume();

  double sum = 0.0;
  const double fac = 4.0 * M_PI * radius * radius * radius / omega;

  // G=0 term
  if ( nspin == 1 )
  {
    sum = fac * cd.rhog[0][0].real() / 3.0;
  }
  else
  {
    // nspin == 2
    if ( ispin == 0 )
      sum = fac * ( cd.rhog[0][0].real() + cd.rhog[1][0].real() ) / 3.0;
    else if ( ispin == 1 )
      sum = fac * cd.rhog[0][0].real() / 3.0;
    else
      sum = fac * cd.rhog[1][0].real() / 3.0;
  }

  // Start sum after G=0 coeff if in first process row
  int igstart = 0;
  if ( ctxt.myrow() == 0 )
  {
    // skip G=0 coefficient
    igstart = 1;
  }

  const int ngloc = cd.rhog[0].size();
  const double *gx = vbasis->gx_ptr(0);
  const double *gy = vbasis->gx_ptr(1);
  const double *gz = vbasis->gx_ptr(2);
  const double *g = vbasis->g_ptr();
  for ( int ig = igstart; ig < ngloc; ig++ )
  {
    // translate origin: compute exp(i G*pos )
    double arg = pos.x*gx[ig] + pos.y*gy[ig] + pos.z*gz[ig];
    double c = cos(arg);
    double s = sin(arg);
    // Re ( c(G) * (c + i*s )) = Re(c(G))*c - Im(c(G)*s

    complex<double> cg;
    if ( nspin == 1 )
    {
      cg = cd.rhog[0][ig];
    }
    else
    {
      // nspin == 2
      if ( ispin == 0 )
        cg = cd.rhog[0][ig] + cd.rhog[1][ig];
      else if ( ispin == 1 )
        cg = cd.rhog[0][ig];
      else
        cg = cd.rhog[1][ig];
    }

    // real part of coefficient of translated function
    double ctrans_re = cg.real() * c - cg.imag() * s;
    // product of norms: g * radius
    double gr = g[ig]*radius;
    // Bessel function j1(z) = sin(z)/z^2 - cos(z)/z
    double j1gr = sin(gr)/(gr*gr) - cos(gr)/gr;
    // factor 2 in next line: G and -G
    sum += 2.0 * fac * ctrans_re * j1gr / gr;
  }

  // accumulate sum across tasks
  ctxt.dsum('c',1,1,&sum,1);

  if ( ui->onpe0() )
  {
    cout << " <partial_charge atom=\"" << atom_name
         << "\" radius=\"" << radius;
    if ( ispin > 0 )
      cout << "\" spin=\"" << ispin;
    cout << "\"> " << setprecision(6) << sum << " </partial_charge>" << endl;
  }

  return 0;
}
