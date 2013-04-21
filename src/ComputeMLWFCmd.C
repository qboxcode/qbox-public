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
// ComputeMLWFCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ComputeMLWFCmd.C,v 1.8 2008-11-14 22:11:30 fgygi Exp $

#include "ComputeMLWFCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;

  // Check that only the k=0 point is used
  if ( wf.nkp()>1 || length(wf.kpoint(0)) != 0.0 )
  {
    if ( ui->onpe0() )
    {
      cout << " ComputeMLWFCmd::action: compute_mlwf can only be used at\n"
           << " the Gamma point (k=0)" << endl;
    }
    return 1;
  }

  if ( ui->onpe0() )
    cout << "<mlwfs>" << endl;

  D3vector edipole_sum;
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    SlaterDet& sd = *(wf.sd(ispin,0));
    MLWFTransform* mlwft = new MLWFTransform(sd);

    mlwft->compute_transform();
    mlwft->apply_transform(sd);

    if ( ui->onpe0() )
    {
      cout << " <mlwfset spin=\"" << ispin
           << "\" size=\"" << sd.nst() << "\">" << endl;
      for ( int i = 0; i < sd.nst(); i++ )
      {
        D3vector ctr = mlwft->center(i);
        double sp = mlwft->spread(i);
        cout.setf(ios::fixed, ios::floatfield);
        cout.setf(ios::right, ios::adjustfield);
        cout << "   <mlwf center=\"" << setprecision(6)
             << setw(12) << ctr.x
             << setw(12) << ctr.y
             << setw(12) << ctr.z
             << " \" spread=\" " << sp << " \"/>"
             << endl;
      }
      D3vector edipole = mlwft->dipole();
      cout << " <electronic_dipole spin=\"" << ispin << "\"> " << edipole
           << " </electronic_dipole>" << endl;
      cout << " </mlwfset>" << endl;
      edipole_sum += edipole;
    }
    delete mlwft;
  }

  if ( ui->onpe0() )
  {
    D3vector idipole = s->atoms.dipole();
    cout << "   <ionic_dipole> " << idipole
         << " </ionic_dipole>" << endl;
    cout << "   <total_dipole> " << idipole + edipole_sum
         << " </total_dipole>" << endl;
    cout << "   <total_dipole_length> " << length(idipole + edipole_sum)
         << " </total_dipole_length>" << endl;
    cout << "</mlwfs>" << endl;
  }
  return 0;
}
