////////////////////////////////////////////////////////////////////////////////
//
// Preconditioner.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Preconditioner.C,v 1.5 2007-10-19 16:24:04 fgygi Exp $

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
  const double ecutpr = s_.ctrl.ecutprec> 0.0 ? s_.ctrl.ecutprec : wf.ecut();

  diag_.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    diag_[ispin].resize(wf.nkp());
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      if ( wf.sd(ispin,ikp) != 0 && wf.sdcontext(ispin,ikp)->active() )
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
}
