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
// PlotCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "PlotCmd.h"
#include "isodate.h"
#include "release.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include "Context.h"
#include "Sample.h"
#include "SampleReader.h"
#include "Basis.h"
#include "EnergyFunctional.h"
#include "FourierTransform.h"
#include "SlaterDet.h"
#include "Matrix.h"
#include "Species.h"
#include "Atom.h"
#include "ChargeDensity.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
int PlotCmd::action(int argc, char **argv)
{
  string usage("  Use: plot filename\n"
               "       plot -density  [-spin {1|2}] filename\n"
               "       plot -vlocal   [-spin {1|2}] filename\n"
               "       plot -wf n [-spin {1|2}] filename\n"
               "       plot -wfs nmin nmax [-spin {1|2}] filename");

  // parse arguments
  // plot filename               : plot atoms in xyz format
  // plot -density filename      : plot atoms and density in cube format
  // plot -vlocal  filename      : plot atoms and vlocal in cube format
  // plot -wf <n> filename       : plot atoms and wf <n> in cube format
  // plot -wf <n1> <n2> filename : plot atoms and wfs <n1> to <n2> in cube fmt
  // spin: 1 = first spin (up), 2 = second spin (down)
  if ( (argc < 2) || (argc > 7) )
  {
    if ( ui->onpe0() )
      cout << usage << endl;
    return 1;
  }

  bool plot_atoms = true;
  bool xyz = true;
  bool plot_density = false;
  bool plot_vlocal = false;
  bool plot_wf = false;
  bool plot_wfs = false;
  int nmin,nmax,nwf;
  // ispin = 0: plot both spins
  // ispin = 1: plot first spin
  // ispin = 2: plot second spin
  int ispin = 0;

  string filename;

  // process arguments
  int iarg = 1;
  while ( iarg < argc )
  {
    if ( !strcmp(argv[iarg],"-density") )
    {
      plot_density = true;
      xyz = false;
    }
    else if ( !strcmp(argv[iarg],"-vlocal") )
    {
      plot_vlocal = true;
      xyz = false;
    }
    else if ( !strcmp(argv[iarg],"-wf") )
    {
      plot_wf = true;
      xyz = false;
      // process argument: n
      iarg++;
      if ( iarg == argc )
      {
        if ( ui->onpe0() )
        cout << usage << endl;
        return 1;
      }
      nmin = atoi(argv[iarg]) - 1;
      nmax = nmin;
      nwf = 1;
      if ( nmin < 0 || nmax >= s->wf.nst() || nmin > nmax )
      {
        if ( ui->onpe0() )
          cout << " nmin or nmax incompatible with nst="
               << s->wf.nst() << endl;
        return 1;
      }
    }
    else if ( !strcmp(argv[iarg],"-wfs") )
    {
      plot_wfs = true;
      xyz = false;
      // process argument: nmin
      iarg++;
      if ( iarg==argc )
      {
        if ( ui->onpe0() )
        cout << usage << endl;
        return 1;
      }
      nmin = atoi(argv[iarg]) - 1;
      // process argument: nmax
      iarg++;
      if ( iarg==argc )
      {
        if ( ui->onpe0() )
        cout << usage << endl;
        return 1;
      }
      nmax = atoi(argv[iarg]) - 1;
      nwf = nmax-nmin+1;
      if ( nmin < 0 || nmax >= s->wf.nst() || nmin > nmax )
      {
        if ( ui->onpe0() )
          cout << " nmin or nmax incompatible with nst="
               << s->wf.nst() << endl;
        return 1;
      }
    }
    else if ( !strcmp(argv[iarg],"-spin") )
    {
      if ( !(plot_density || plot_wf || plot_wfs || plot_vlocal) )
      {
        if ( ui->onpe0() )
          cout << usage << endl;
        return 1;
      }
      if ( s->wf.nspin() != 2 )
      {
        if ( ui->onpe0() )
          cout << "nspin = 1, cannot select spin" << endl;
        return 1;
      }
      // process argument: ispin
      iarg++;
      if ( iarg==argc )
      {
        if ( ui->onpe0() )
          cout << usage << endl;
        return 1;
      }
      ispin = atoi(argv[iarg]);
      if ( ispin < 1 || ispin > 2 )
      {
        if ( ui->onpe0() )
          cout << " spin must be 1 or 2" <<  endl;
        return 1;
      }
    }
    else
    {
      // argument must be the file name
      filename = argv[iarg];
    }

    iarg++;

  } // while iarg

  // Must specify spin if plotting wave functions when nspin==2
  if ( s->wf.nspin()==2 && (plot_vlocal||plot_wf||plot_wfs) && ispin==0 )
  {
    if ( ui->onpe0() )
      cout << " must use -spin if nspin==2" <<  endl;
    return 1;
  }

  ofstream os;
  vector<double> tmpr;
  int np0 = 0;
  int np1 = 0;
  int np2 = 0;

  const Context& ctxt = s->wf.sd_context();
  if ( plot_density )
  {
    ChargeDensity cd(s->wf);
    cd.update_density();
    tmpr.resize(cd.vft()->np012());
    np0 = cd.vft()->np0();
    np1 = cd.vft()->np1();
    np2 = cd.vft()->np2();
    if ( s->wf.nspin() == 1 )
    {
      for ( int i = 0; i < cd.vft()->np012loc(); i++ )
        tmpr[i] = cd.rhor[0][i];
    }
    else
    {
      if ( ispin == 0 )
      {
        // plot both spins
        for ( int i = 0; i < cd.vft()->np012loc(); i++ )
          tmpr[i] = cd.rhor[0][i] + cd.rhor[1][i];
      }
      else
      {
        // plot one spin only
        // ispin==1 or ispin==2
        for ( int i = 0; i < cd.vft()->np012loc(); i++ )
          tmpr[i] = cd.rhor[ispin-1][i];
      }
    }

    // send blocks of tmpr to pe0
    // send from first context column only
    for ( int i = 0; i < ctxt.nprow(); i++ )
    {
      bool iamsending = ctxt.mycol() == 0 && i == ctxt.myrow();

      // send size of tmpr block
      int size=-1;
      if ( ctxt.onpe0() )
      {
        if ( iamsending )
        {
          // sending to self, size not needed
        }
        else
          ctxt.irecv(1,1,&size,1,i,0);
      }
      else
      {
        if ( iamsending )
        {
          size = cd.vft()->np012loc();
          ctxt.isend(1,1,&size,1,0,0);
        }
      }

      // send tmpr block
      if ( ctxt.onpe0() )
      {
        if ( iamsending )
        {
          // do nothing, data is already in place
        }
        else
        {
          int istart = cd.vft()->np0() * cd.vft()->np1() *
                       cd.vft()->np2_first(i);
          ctxt.drecv(size,1,&tmpr[istart],1,i,0);
        }
      }
      else
      {
        if ( iamsending )
        {
          ctxt.dsend(size,1,&tmpr[0],1,0,0);
        }
      }
    }
  } // plot_density
  else if ( plot_vlocal )
  {
    ChargeDensity cd(s->wf);
    EnergyFunctional ef(*s,cd);
    cd.update_density();
    cd.update_rhor();
    bool compute_stress = false;
    ef.update_vhxc(compute_stress);

    tmpr.resize(cd.vft()->np012());
    np0 = cd.vft()->np0();
    np1 = cd.vft()->np1();
    np2 = cd.vft()->np2();
    if ( s->wf.nspin() == 1 )
    {
      for ( int i = 0; i < cd.vft()->np012loc(); i++ )
        tmpr[i] = ef.v_r[0][i];
    }
    else
    {
      if ( ispin == 0 )
      {
        // should not get here
        // ispin must be set if nspin==2
        if ( ui->onpe0() )
          cout << usage << endl;
        return 1;
      }
      else
      {
        // plot one spin only
        // ispin==1 or ispin==2
        for ( int i = 0; i < cd.vft()->np012loc(); i++ )
          tmpr[i] = ef.v_r[ispin-1][i];
      }
    }

    // send blocks of tmpr to pe0
    // send from first context column only
    for ( int i = 0; i < ctxt.nprow(); i++ )
    {
      bool iamsending = ctxt.mycol() == 0 && i == ctxt.myrow();

      // send size of tmpr block
      int size=-1;
      if ( ctxt.onpe0() )
      {
        if ( iamsending )
        {
          // sending to self, size not needed
        }
        else
          ctxt.irecv(1,1,&size,1,i,0);
      }
      else
      {
        if ( iamsending )
        {
          size = cd.vft()->np012loc();
          ctxt.isend(1,1,&size,1,0,0);
        }
      }

      // send tmpr block
      if ( ctxt.onpe0() )
      {
        if ( iamsending )
        {
          // do nothing, data is already in place
        }
        else
        {
          int istart = cd.vft()->np0() * cd.vft()->np1() *
                       cd.vft()->np2_first(i);
          ctxt.drecv(size,1,&tmpr[istart],1,i,0);
        }
      }
      else
      {
        if ( iamsending )
        {
          ctxt.dsend(size,1,&tmpr[0],1,0,0);
        }
      }
    }
  } // plot_vlocal
  else if ( plot_wf || plot_wfs )
  {
    // compute wf or wf squared and store in tmpr
    if ( ctxt.onpe0() )
    {
      ctxt.ibcast_send(1,1,&nwf,1);
      ctxt.ibcast_send(1,1,&nmin,1);
      ctxt.ibcast_send(1,1,&nmax,1);
    }
    else
    {
      ctxt.ibcast_recv(1,1,&nwf,1,0,0);
      ctxt.ibcast_recv(1,1,&nmin,1,0,0);
      ctxt.ibcast_recv(1,1,&nmax,1,0,0);
    }

    if ( nwf > 0 && s->wf.nst() == 0 )
    {
      cout << " no states in sample" << endl;
      return 1;
    }

    int isp_min, isp_max;
    SlaterDet *sdp;
    if ( ispin == 0 )
    {
      // -spin was not specified:
      if ( s->wf.nspin() == 1 )
      {
        isp_min = 0; isp_max = 0;
      }
      else
      {
        isp_min = 0; isp_max = 1;
      }
    }
    else
    {
      isp_min = ispin-1; isp_max = ispin-1;
    }

    const Basis& basis = s->wf.sd(0,0)->basis();
    np0 = basis.np(0);
    np1 = basis.np(1);
    np2 = basis.np(2);
    FourierTransform ft(basis,np0,np1,np2);

    vector<complex<double> > wftmp(ft.np012loc());
    vector<double> wftmpr(ft.np012());
    tmpr.resize(ft.np012());
    tmpr.assign(ft.np012(),0.0);

    for ( int isp = isp_min; isp <= isp_max; isp++ )
    {
      sdp = s->wf.sd(isp,0);
      const ComplexMatrix& c = sdp->c();

      for ( int n = nmin; n <= nmax; n++ )
      {
        if ( n >= s->wf.nst(isp) )
        {
          if ( ui->onpe0() )
            cout << "invalid wave function index: " << n+1
                 << " > nst(ispin)" << endl;
          return 1;
        }

        // compute real-space wavefunction

        // transform wf on ctxt.mycol() hosting state n
        if ( c.pc(n) == c.context().mycol() )
        {
          //os << " state " << n << " is stored on column "
          //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
          int nloc = c.y(n); // local index
          ft.backward(c.cvalptr(c.mloc()*nloc),&wftmp[0]);

          double *a = (double*) &wftmp[0];
          if ( basis.real() )
          {
            // real function: plot wf
            for ( int i = 0; i < ft.np012loc(); i++ )
              wftmpr[i] = a[2*i];
          }
          else
          {
            // complex function: plot modulus
            for ( int i = 0; i < ft.np012loc(); i++ )
              wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
          }
        }

        // send blocks of wftmpr to pe0
        for ( int i = 0; i < c.context().nprow(); i++ )
        {
          bool iamsending = c.pc(n) == c.context().mycol() &&
                            i == c.context().myrow();

          // send size of wftmpr block
          int size=-1;
          if ( c.context().onpe0() )
          {
            if ( iamsending )
            {
              // sending to self, size not needed
            }
            else
              c.context().irecv(1,1,&size,1,i,c.pc(n));
          }
          else
          {
            if ( iamsending )
            {
              size = ft.np012loc();
              c.context().isend(1,1,&size,1,0,0);
            }
          }

          // send wftmpr block
          if ( c.context().onpe0() )
          {
            if ( iamsending )
            {
              // do nothing, data is already in place
            }
            else
            {
              int istart = ft.np0() * ft.np1() * ft.np2_first(i);
              c.context().drecv(size,1,&wftmpr[istart],1,i,c.pc(n));
            }
          }
          else
          {
            if ( iamsending )
            {
              c.context().dsend(size,1,&wftmpr[0],1,0,0);
            }
          }
        }

        // process the data on task 0
        if ( c.context().onpe0() )
        {
          // wftmpr is now complete on task 0
          if ( plot_wfs )
          {
            // multiple wfs, accumulate square
            for ( int i = 0; i < ft.np012(); i++ )
            {
              tmpr[i] += wftmpr[i]*wftmpr[i];
            }
          }
          else
          {
            // plot individual wf
            for ( int i = 0; i < ft.np012(); i++ )
            {
              tmpr[i] = wftmpr[i];
            }
          }
        }
      } // for n
    } // for isp
  } // if plot_wf || plot_wfs

  // tmpr now contains the function to plot on task 0

  if ( ctxt.onpe0() )
    os.open(filename.c_str());

  if ( plot_atoms )
  {
    if ( ctxt.onpe0() )
    {
      if ( xyz )
      {
        const double a0 = 0.529177;
        int natoms = s->atoms.size();
        os << natoms << endl;
        os << "Created " << isodate() << " by qbox-" << release() << endl;
        const int nsp = s->atoms.nsp();
        for ( int is = 0; is < nsp; is++ )
        {
          Species* sp = s->atoms.species_list[is];
          string symbol = sp->symbol();
          const int na = s->atoms.na(is);
          for ( int ia = 0; ia < na; ia++ )
          {
            Atom *ap = s->atoms.atom_list[is][ia];
            os << setprecision(5);
            os << symbol << " " << a0*ap->position() << endl;
          }
        }
      }
      else
      {
        // write header and atoms
        os << "Created " << isodate() << " by qbox-" << release() << endl;
        os << endl;

        int natoms = s->atoms.size();
        D3vector a0 = s->atoms.cell().a(0);
        D3vector a1 = s->atoms.cell().a(1);
        D3vector a2 = s->atoms.cell().a(2);
        os << natoms << " " << -0.5*(a0+a1+a2) << endl;

        // write unit cell
        os << np0 << " " << a0/np0 << endl;
        os << np1 << " " << a1/np1 << endl;
        os << np2 << " " << a2/np2 << endl;
        const int nsp = s->atoms.nsp();
        for ( int is = 0; is < nsp; is++ )
        {
          Species* sp = s->atoms.species_list[is];
          const int z = sp->atomic_number();
          const int na = s->atoms.na(is);
          for ( int ia = 0; ia < na; ia++ )
          {
            Atom *ap = s->atoms.atom_list[is][ia];
            os << setprecision(5);
            os << z << " " << ((double) z) << " " << ap->position() << endl;
          }
        }
      }
    }
  } // if plot_atoms

  if ( plot_density || plot_vlocal || plot_wf || plot_wfs )
  {
    // process the function in tmpr
    if ( ctxt.onpe0() )
    {
      os.setf(ios::scientific,ios::floatfield);
      os << setprecision(5);
      for ( int i = 0; i < np0; i++ )
      {
        const int ip = (i + np0/2 ) % np0;
        for ( int j = 0; j < np1; j++ )
        {
          const int jp = (j + np1/2 ) % np1;
          for ( int k = 0; k < np2; k++ )
          {
            const int kp = (k + np2/2 ) % np2;
            os << setw(13) << tmpr[ip+np0*(jp+np1*kp)];
            if ( ( k % 6 ) == 5 )
              os << endl;
          }
          if ( ( np2 % 6 ) != 0 )
            os << endl;
        }
      }
    }
  } // if plot_density || ...

  os.close();

  return 0;
}
