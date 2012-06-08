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
// Preconditioner.C
//
////////////////////////////////////////////////////////////////////////////////

#include "Preconditioner.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Preconditioner::Preconditioner(const Wavefunction& wf, double ecutprec)
 : ecutprec_(ecutprec)
{
  kpg2_.resize(wf.nspin());
  ekin_.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    kpg2_[ispin].resize(wf.nkp());
    ekin_[ispin].resize(wf.nkp());
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      const SlaterDet& sd = *(wf.sd(ispin,ikp));
      const Basis& wfbasis = sd.basis();
      kpg2_[ispin][ikp] = wfbasis.kpg2_ptr();
      ekin_[ispin][ikp].resize(sd.nstloc());
    }
  }
  update(wf);
}

////////////////////////////////////////////////////////////////////////////////
double Preconditioner::diag(int ispin, int ikp, int n, int ig) const
{
  if ( ecutprec_ == 0.0 )
  {
    const double ekin_n = ekin_[ispin][ikp][n];
    if ( ekin_n == 0.0 )
      return 0.0;
    const double q2 = kpg2_[ispin][ikp][ig];
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
  }
  else
  {
    // next lines: inclusion of fstress if confinement potential is used
    // const valarray<double>& fstress = ef_.confpot(ikp)->fstress();
    // double e = 0.5 * ( kpg2_ptr[ig] + fstress[ig] );
    // diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;

    const double e = 0.5 * kpg2_[ispin][ikp][ig];
    return ( e < ecutprec_ ) ? 0.5 / ecutprec_ : 0.5 / e;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Preconditioner::update(const Wavefunction& wf)
{
  // update the kinetic energy ekin_[ispin][ikp][n] of states in wf
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      const SlaterDet& sd = *(wf.sd(ispin,ikp));
      const Basis& wfbasis = sd.basis();
      // factor fac in next lines: 2.0 for G and -G (if basis is real) and
      // 0.5 from 1/(2m)
      const double fac = wfbasis.real() ? 1.0 : 0.5;

      const ComplexMatrix& c = sd.c();
      const Context& sdctxt = sd.context();

      const int ngwloc = wfbasis.localsize();
      const complex<double>* p = c.cvalptr();
      const int mloc = c.mloc();
      const int nloc = c.nloc();

      valarray<double> buf(2*nloc);
      const double* pkpg2 = kpg2_[ispin][ikp];
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
      for ( int n = 0; n < nloc; n++ )
        ekin_[ispin][ikp][n] = buf[n] > 0.0 ? buf[n+nloc]/buf[n] : 0;

#ifdef DEBUG
      if ( sdctxt.onpe0() )
      {
        for ( int n = 0; n < nloc; n++ )
          cout << "Preconditioner::update ekin[" << n << "] = "
               << ekin_[ispin][ikp][n] << endl;
      }
#endif
    }
  }
}
