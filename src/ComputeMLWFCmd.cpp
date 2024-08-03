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
// ComputeMLWFCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "ComputeMLWFCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
#include "cout0.h"
#include <stdexcept>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;

  const bool onpe0 = MPIdata::onpe0();

  // Check that only the k=0 point is used
  if ( wf.nkp()>1 || length(wf.kpoint(0)) != 0.0 )
    throw runtime_error("ComputeMLWFCmd: can only be used with 1 k-point");

  // Check if nempty > 0
  if ( wf.nempty() > 0 )
    throw runtime_error("ComputeMLWFCmd: cannot run with nempty > 0");

  if ( onpe0 )
    cout << "<mlwfs>" << endl;

  D3vector edipole_sum;
  for ( int ispin = 0; ispin < wf.nspin(); ++ispin )
  {
    const int isp_loc = wf.isp_local(ispin);
    const int ikp = 0;
    const int ikp_loc = wf.ikp_local(ikp);
    ostringstream ostr;
    int isrc = -1;
    if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
    {
      SlaterDet& sd = *(wf.sd(isp_loc,ikp_loc));
      MLWFTransform* mlwft = new MLWFTransform(sd);

      mlwft->update();
      mlwft->compute_transform();
      mlwft->apply_transform(sd);

      if ( MPIdata::sd_rank() == 0 )
      {
        ostr.str("");
        isrc = MPIdata::rank();
        ostr << " <mlwfset spin=\"" << ispin
             << "\" size=\"" << sd.nst() << "\">" << endl;
        edipole_sum += mlwft->dipole();
        ostr << *mlwft;
        ostr << " </mlwfset>" << endl;
      } // sd_rank() == 0
      delete mlwft;
    }
    cout0(ostr.str(),isrc);
    MPI_Barrier(MPIdata::comm());
  } // ispin
  if ( onpe0 )
    cout << "</mlwfs>" << endl;

  if ( MPIdata::sd_rank() == 0 )
  {
    D3vector d3tsum;
    MPI_Reduce(&edipole_sum[0],&d3tsum[0],3,
               MPI_DOUBLE,MPI_SUM,0,MPIdata::sp_comm());
    edipole_sum = d3tsum;
  }

  if ( onpe0 )
  {
    D3vector idipole = s->atoms.dipole();
    cout << setprecision(8) << fixed << right;
    cout << "<dipole>" << endl;
    cout << " <dipole_ion>   "
         << setw(14) << idipole.x
         << setw(14) << idipole.y
         << setw(14) << idipole.z
         << " </dipole_ion>" << endl;
    cout << " <dipole_el>    "
         << setw(14) << edipole_sum.x
         << setw(14) << edipole_sum.y
         << setw(14) << edipole_sum.z
         << " </dipole_el>" << endl;
    D3vector dipole_total = idipole + edipole_sum;
    cout << " <dipole_total> "
         << setw(14) << dipole_total.x
         << setw(14) << dipole_total.y
         << setw(14) << dipole_total.z
         << " </dipole_total>" << endl;
    cout << "</dipole>" << endl;
  }
  return 0;
}
