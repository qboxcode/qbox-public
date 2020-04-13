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
// SpectrumCmd.C:
//
////////////////////////////////////////////////////////////////////////////////

#include "SpectrumCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
int SpectrumCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;

  // Check that only the k=0 point is used
  if ( wf.nkp()>1 || length(wf.kpoint(0)) != 0.0 )
  {
    if ( ui->onpe0() )
    {
      cout << " SpectrumCmd::action: spectrum can only be used at\n"
           << " the Gamma point (k=0)" << endl;
    }
    return 1;
  }

  const UnitCell& cell = wf.cell();
  if ( ui->onpe0() )
    cout << "<dipole_matrix_elements>" << endl;

  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    SlaterDet& sd = *(wf.sd(ispin,0));
    MLWFTransform* mlwft = new MLWFTransform(sd);

    mlwft->update();

    if ( ui->onpe0() )
    {
      const int n = sd.nst();
      const double *c0 = mlwft->a(0)->cvalptr();
      const double *s0 = mlwft->a(1)->cvalptr();
      const double *c1 = mlwft->a(2)->cvalptr();
      const double *s1 = mlwft->a(3)->cvalptr();
      const double *c2 = mlwft->a(4)->cvalptr();
      const double *s2 = mlwft->a(5)->cvalptr();
      const double fac = 0.5 / M_PI;

      cout.setf(ios::fixed, ios::floatfield);
      cout.setf(ios::right, ios::adjustfield);

      for ( int i = 0; i < n; i++ )
      {
        for ( int j = i; j < n; j++ )
        {
          double c[3] = { fac*c0[i+n*j], fac*c1[i+n*j], fac*c2[i+n*j] };
          double s[3] = { fac*s0[i+n*j], fac*s1[i+n*j], fac*s2[i+n*j] };

          // cc, ss: matrix elements in cartesian coordinates
          double cc[3], ss[3];
          cell.vecmult3x3(cell.amat(),c,cc);
          cell.vecmult3x3(cell.amat(),s,ss);

          cout << " <dipole i=\"" << i+1 << "\" j=\"" << j+1 << "\">" << endl;
          cout << setw(12) << setprecision(6) << cc[0] << " "
               << setw(12) << setprecision(6) << ss[0] << endl;
          cout << setw(12) << setprecision(6) << cc[1] << " "
               << setw(12) << setprecision(6) << ss[1] << endl;
          cout << setw(12) << setprecision(6) << cc[2] << " "
               << setw(12) << setprecision(6) << ss[2] << endl;
          cout << cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2]+
                  ss[0]*ss[0]+ss[1]*ss[1]+ss[2]*ss[2] << endl;
          cout << " </dipole>" << endl;
        }
      }
    }
    delete mlwft;
  }

  if ( ui->onpe0() )
    cout << "</dipole_matrix_elements>" << endl;
  return 0;
}
