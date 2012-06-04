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
// $Id: Preconditioner.C,v 1.8 2008-09-08 15:56:19 fgygi Exp $

#include "Preconditioner.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ConfinementPotential.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
Preconditioner::Preconditioner(const Sample& s, const EnergyFunctional& ef) :
  s_(s), ef_(ef)
{
  update();
}

////////////////////////////////////////////////////////////////////////////////
void Preconditioner::update(void)
{
  // reinitialize preconditioner
  bool use_confinement = s_.ctrl.ecuts > 0.0;
  const Wavefunction& wf = s_.wf;
  // If ecutprec is zero, use ecut
  const double ecutpr = s_.ctrl.ecutprec > 0.0 ? s_.ctrl.ecutprec : wf.ecut();

  diag_.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    diag_[ispin].resize(wf.nkp());
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      // Only resize and initialize diag_ if ikp is active on this task
      const Basis& basis = wf.sd(ispin,ikp)->basis();
      const int ngwloc = basis.localsize();
      diag_[ispin][ikp].resize(ngwloc);
      const double *kpg2_ptr = basis.kpg2_ptr();

      if ( use_confinement )
      {
        const valarray<double>& fstress = ef_.confpot(ikp)->fstress();
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          double e = 0.5 * ( kpg2_ptr[ig] + fstress[ig] );
          diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
        }
      }
      else
      {
        for ( int ig = 0; ig < ngwloc; ig++ )
        {
          double e = 0.5 * kpg2_ptr[ig];
          diag_[ispin][ikp][ig] = ( e < ecutpr ) ? 0.5 / ecutpr : 0.5 / e;
        }
      }
    }
  }
}
