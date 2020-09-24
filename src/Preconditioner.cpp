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
// Preconditioner.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Preconditioner.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "EnergyFunctional.h"
#include "ConfinementPotential.h"
#include "SlaterDet.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Preconditioner::Preconditioner(const Wavefunction& wf, EnergyFunctional& ef,
  double ecutprec) : ef_(ef), ecutprec_(ecutprec)
{
  kpg2_.resize(wf.nsp_loc());
  ekin_.resize(wf.nsp_loc());
  for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
  {
    kpg2_[isp_loc].resize(wf.nkp_loc());
    ekin_[isp_loc].resize(wf.nkp_loc());
    for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
    {
      const SlaterDet& sd = *(wf.sd(isp_loc,ikp_loc));
      const Basis& wfbasis = sd.basis();
      kpg2_[isp_loc][ikp_loc] = wfbasis.kpg2_ptr();
      ekin_[isp_loc][ikp_loc].resize(sd.nstloc());
    }
  }
  update(wf);
}

////////////////////////////////////////////////////////////////////////////////
double Preconditioner::diag(int isp_loc, int ikp_loc, int n, int ig) const
{
  const valarray<double>& fstress = ef_.confpot(ikp_loc)->fstress();
  if ( ecutprec_ == 0.0 )
  {
    double ekin_n = ekin_[isp_loc][ikp_loc][n];
    // if ekin_n == 0 (occurs for first wf, G=0, when starting without
    // randomizing wfs) replace ekin_n by fixed value 1.0
    if ( ekin_n == 0.0 )
      ekin_n = 1.0;
#if 0
    const double q2 = kpg2_[isp_loc][ikp_loc][ig] + fstress[ig];
#if 0
    // Payne Teter Allan adaptive preconditioner
    const double x = 0.5 * q2 / ekin_n;
    const double num = 27.0 + x * ( 18.0 + x * ( 12 + x * 8 ));
    const double den = num + 16.0 * x*x*x*x;
    return (num/den);
#else
    // modified Payne Teter Allan adaptive preconditioner (Kresse)
    const double x = 0.5 * q2 / ( 1.5 * ekin_n );
    const double num = 27.0 + x * ( 18.0 + x * ( 12 + x * 8 ));
    const double den = num + 16.0 * x*x*x*x;
    return ( 2.0 / (1.5 *  ekin_n )) * (num/den);
#endif
#else
    // basic adaptive preconditioner: use ekin_n for the value of ecutprec
    double e = 0.5 * ( kpg2_[isp_loc][ikp_loc][ig] + fstress[ig] );
    return ( e < ekin_n ) ? 0.5 / ekin_n : 0.5 / e;
#endif
  }
  else
  {
    double e = 0.5 * ( kpg2_[isp_loc][ikp_loc][ig] + fstress[ig] );
    return ( e < ecutprec_ ) ? 0.5 / ecutprec_ : 0.5 / e;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Preconditioner::update(const Wavefunction& wf)
{
  // update the kinetic energy ekin_[isp_loc][ikp_loc][n] of states in wf
  for ( int isp_loc = 0; isp_loc < wf.nsp_loc(); ++isp_loc )
  {
    for ( int ikp_loc = 0; ikp_loc < wf.nkp_loc(); ++ikp_loc )
    {
      const SlaterDet& sd = *(wf.sd(isp_loc,ikp_loc));
      const Basis& wfbasis = sd.basis();
      // factor fac in next lines: 2.0 for G and -G (if basis is real)
      const double fac = wfbasis.real() ? 2.0 : 1.0;

      const ComplexMatrix& c = sd.c();
      const Context& sdctxt = sd.context();

      const int ngwloc = wfbasis.localsize();
      const complex<double>* p = c.cvalptr();
      const int mloc = c.mloc();
      const int nloc = c.nloc();

      valarray<double> buf(2*nloc);
      const double* pkpg2 = kpg2_[isp_loc][ikp_loc];
      for ( int n = 0; n < nloc; n++ )
      {
        double sum_norm = 0.0;
        double sum_ekin = 0.0;
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          const double psi2 = norm(p[ig+n*mloc]);
          sum_norm += psi2;
          sum_ekin += psi2 * pkpg2[ig];
        }
        // correct for double counting of G=0 in norm
        if ( sdctxt.myrow() == 0 )
          sum_norm -= 0.5 * norm(p[n*mloc]);
        // store norm in buf[n] and ekin in buf[n+nloc]
        buf[n] = fac * sum_norm;
        buf[n+nloc] = fac * sum_ekin;
      }
      sdctxt.dsum('C',2*nloc,1,&buf[0],2*nloc);
      // factor 0.5 in next line: 1/(2m)
      for ( int n = 0; n < nloc; n++ )
        ekin_[isp_loc][ikp_loc][n] =
          0.5*buf[n] > 0.0 ? 0.5*buf[n+nloc]/buf[n] : 0;

#ifdef DEBUG
      if ( sdctxt.onpe0() )
      {
        for ( int n = 0; n < nloc; n++ )
          cout << "Preconditioner::update ekin[" << n << "] = "
               << ekin_[isp_loc][ikp_loc][n] << endl;
      }
#endif
    }
  }
}
