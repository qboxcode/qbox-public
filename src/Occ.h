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
// Occ.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef OCC_H
#define OCC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Occ : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "occ"; };

  int set ( int argc, char **argv )
  {
    if ( (argc!=3) && (argc!=4) )
    {
      if ( ui->onpe0() )
      {
        cout << " use: set occ [ispin] n f" << endl;
        cout << " ispin = {1|2}" << endl;
      }
      return 1;
    }

    if ( s->wf.nkp() != 1 )
    {
      if ( ui->onpe0() )
        cout << " set occ: not implemented for multiple k-points" << endl;
      return 1;
    }

    if ( argc == 3 )
    {
      // set occ n f
      if ( s->wf.nspin() != 1 )
      {
        if ( ui->onpe0() )
          cout << "nspin must be 1" << endl;
        return 1;
      }

      const int n = atoi(argv[1]);
      if ( (n < 1) || n > s->wf.nst() )
      {
        if ( ui->onpe0() )
          cout << " n must be in [1," << s->wf.nst() << "]" << endl;
        return 1;
      }

      const double f = atof(argv[2]);
      if ( (f < 0.0) || (f > 2.0) )
      {
        if ( ui->onpe0() )
          cout << " f must be in [0,2]" << endl;
        return 1;
      }

      Wavefunction& wf = s->wf;
      SlaterDet& sd = *wf.sd(0,0);
      // n-1 in next line: states are numbered starting from 1 in
      // the set occ command
      sd.set_occ(n-1,f);
      // recompute total electronic charge
      double sum = 0.0;
      for ( int i = 0; i < sd.nst(); i++ )
        sum += sd.occ(i);
      if ( ui->onpe0() )
        cout << " total electronic charge: " << sum << endl;
      // adjust total number of electrons
      //wf.set_nel((int)sum);
    }

    if ( argc == 4 )
    {
      // set occ ispin n f
      if ( s->wf.nspin() != 2 )
      {
        if ( ui->onpe0() )
          cout << "nspin must be 2 when using ispin" << endl;
        return 1;
      }

      const int ispin  = atoi(argv[1]);
      if ( (ispin < 1) || (ispin > 2) )
      {
        if ( ui->onpe0() )
          cout << " ispin must be 1 or 2" << endl;
        return 1;
      }

      const int n = atoi(argv[2]);
      if ( (n < 1) || n > s->wf.nst(ispin-1) )
      {
        if ( ui->onpe0() )
          cout << " n must be in [1," << s->wf.nst(ispin-1) << "]" << endl;
        return 1;
      }

      const double f = atof(argv[3]);
      if ( (f < 0.0) || (f > 1.0) )
      {
        if ( ui->onpe0() )
          cout << " f must be in [0,1] when using ispin" << endl;
        return 1;
      }

      Wavefunction& wf = s->wf;
      // ispin-1 in next line: spins are numbered starting from 1 in
      // the set occ command
      SlaterDet& sd = *wf.sd(ispin-1,0);
      // n-1 in next line: states are numbered starting from 1 in
      // the set occ command
      sd.set_occ(n-1,f);
      // recompute total electronic charge
      double sum = 0.0;
      for ( int isp = 0; isp < wf.nspin(); isp++ )
        for ( int i = 0; i < wf.nst(isp); i++ )
          sum += wf.sd(isp,0)->occ(i);
      if ( ui->onpe0() )
        cout << " total electronic charge: " << sum << endl;
      // adjust total number of electrons
      //wf.set_nel((int)sum);
    }
    return 0;
  }

  string print (void) const
  {
    ostringstream st;
    st.setf(ios::right,ios::adjustfield);
    st.setf(ios::fixed,ios::floatfield);

    const Wavefunction& wf = s->wf;
    st << " occupation numbers" << endl;
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      const SlaterDet& sd = *wf.sd(ispin,0);
      const int nst = sd.nst();
      for ( int n = 0; n < nst; n++ )
      {
        st << setw(7) << setprecision(4) << sd.occ(n);
        if ( ( n%10 ) == 9 ) st << endl;
      }
      if ( nst % 10 != 0 )
        st << endl;
      if ( (wf.nspin() == 2) && (ispin == 0) )
        st << endl;
    }
    st << " total electronic charge: " << wf.nel();
    return st.str();
 }

  Occ(Sample *sample) : s(sample) {};
};
#endif
