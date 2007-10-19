////////////////////////////////////////////////////////////////////////////////
//
// ComputeMLWFCmd.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ComputeMLWFCmd.C,v 1.5 2007-10-19 16:24:04 fgygi Exp $

#include "ComputeMLWFCmd.h"
#include<iostream>
#include "Context.h"
#include "SlaterDet.h"
using namespace std;

int ComputeMLWFCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;
  SlaterDet& sd = *(wf.sd(0,0));

  mlwft = new MLWFTransform(sd);

  mlwft->compute_transform();
  mlwft->apply_transform(sd);

  if ( ui->onpe0() )
  {
    cout << " <mlwf_set size=\"" << sd.nst() << "\">" << endl;
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
    cout << " </mlwf_set>" << endl;
    D3vector edipole = mlwft->dipole();
    cout << " <electronic_dipole> " << edipole
         << " </electronic_dipole>" << endl;
    D3vector idipole = s->atoms.dipole();
    cout << " <ionic_dipole> " << idipole
         << " </ionic_dipole>" << endl;
    cout << " <total_dipole> " << idipole + edipole
         << " </total_dipole>" << endl;
    cout << " <total_dipole_length> " << length(idipole + edipole)
         << " </total_dipole_length>" << endl;
  }
  delete mlwft;
  return 0;
}
