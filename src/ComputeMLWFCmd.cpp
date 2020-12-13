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
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;

  const bool onpe0 = MPIdata::onpe0();

  // Check that only the k=0 point is used
  if ( wf.nkp()>1 || length(wf.kpoint(0)) != 0.0 )
  {
    if ( onpe0 )
    {
      cout << " ComputeMLWFCmd::action: compute_mlwf can only be used at\n"
           << " the Gamma point (k=0)" << endl;
    }
    return 1;
  }

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
        double total_spread[6];
        for ( int j = 0; j < 6; j++ )
           total_spread[j] = 0.0;
        for ( int i = 0; i < sd.nst(); i++ )
        {
          D3vector ctr = mlwft->center(i);
          double spi[6];
          for (int j=0; j<3; j++)
          {
            spi[j] = mlwft->spread2(i,j);
            total_spread[j] += spi[j];
          }

          ostr.setf(ios::fixed, ios::floatfield);
          ostr.setf(ios::right, ios::adjustfield);
          ostr << "   <mlwf center=\"" << setprecision(6)
               << setw(12) << ctr.x
               << setw(12) << ctr.y
               << setw(12) << ctr.z
               << " \" spread=\" "
               << setw(12) << spi[0]
               << setw(12) << spi[1]
               << setw(12) << spi[2] << " \"/>"
               << endl;
        }

        ostr << " <total_spread> ";
        for ( int j = 0; j < 3; j++ )
          ostr << setprecision(6) << setw(15) << total_spread[j];
        ostr << " </total_spread>" << endl;
        D3vector edipole = mlwft->dipole();
        ostr << " <electronic_dipole spin=\"" << ispin
             << "\"> " << edipole << " </electronic_dipole>" << endl;
        edipole_sum += edipole;
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
    cout << setprecision(6);
    cout << "   <ionic_dipole> " << idipole
         << " </ionic_dipole>" << endl;
    cout << "   <total_dipole> " << idipole + edipole_sum
         << " </total_dipole>" << endl;
    cout << "   <total_dipole_length> " << length(idipole + edipole_sum)
         << " </total_dipole_length>" << endl;
  }
  return 0;
}
