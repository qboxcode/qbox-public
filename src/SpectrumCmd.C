////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2020 The Regents of the University of California
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
#include <iostream>
#include <fstream>
#include <cassert>
#include "Context.h"
#include "ChargeDensity.h"
#include "EnergyFunctional.h"
#include "MLWFTransform.h"
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
      cout << " SpectrumCmd::action: spectrum only implemented at\n"
           << " the Gamma point (k=0)" << endl;
    }
    return 1;
  }

  if ( !( argc == 2 || argc == 3 || argc == 5 ) )
  {
    if ( ui->onpe0() )
    {
      cout << " use: spectrum [emin emax] [width] filename" << endl;
    }
    return 1;
  }

  // Compute eigenvalues using the current wave function wf
  Wavefunction dwf(wf);
  ChargeDensity cd(wf);
  cd.update_density();
  EnergyFunctional ef(*s,cd);
  const bool compute_stress = false;
  ef.update_vhxc(compute_stress);
  const bool compute_forces = false;
  const bool compute_hpsi = true;
  valarray<double> sigma_eks;
  vector<vector<double> > fion;
  ef.energy(compute_hpsi,dwf,compute_forces,fion,compute_stress,sigma_eks);
  const bool compute_eigvec = true;
  wf.diag(dwf,compute_eigvec);

  if ( ui->onpe0() )
  {
    cout << "<eigenset>" << endl;
    // print eigenvalues
    for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
    {
      for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
      {
        const int nst = wf.sd(ispin,ikp)->nst();
        const double eVolt = 2.0 * 13.6058;
        cout <<    "  <eigenvalues spin=\"" << ispin
             << "\" kpoint=\""
             << setprecision(8)
             << wf.sd(ispin,ikp)->kpoint()
             << "\" weight=\""
             << setprecision(8)
             << wf.weight(ikp)
             << "\" n=\"" << nst << "\">" << endl;
        for ( int i = 0; i < nst; i++ )
        {
          cout << setw(12) << setprecision(5)
               << wf.sd(ispin,ikp)->eig(i)*eVolt;
          if ( i%5 == 4 ) cout << endl;
        }
        if ( nst%5 != 0 ) cout << endl;
        cout << "  </eigenvalues>" << endl;
      }
    }
    cout << "</eigenset>" << endl;
  }

  const double eVolt = 2.0 * 13.6058;
  const UnitCell& cell = wf.cell();

  // emin, emax: bounds of plot in eV
  // de: energy spacing of plot values in eV (fixed at 0.01)
  // width: gaussian width of convolution in eV (default 0.05)
  const double de = 0.01;
  double emin = 0.0, emax = 0.0, width = 0.05, erange = 0.0;
  const char *spfilename;

  // spectrum file
  if ( argc == 2 )
  {
    spfilename = argv[1];
  }

  // spectrum width file
  if ( argc == 3 )
  {
    width = atof(argv[1]);
    spfilename = argv[2];
  }

  // spectrum emin emax de width file
  if ( argc == 5 )
  {
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    width = atof(argv[3]);
    spfilename = argv[4];
    erange = emax - emin + 3 * width;
    if ( emax <= emin )
    {
      if ( ui->onpe0() )
      {
        cout << " SpectrumCmd::action: emax must be larger than emin" << endl;
      }
      return 1;
    }
  }

  const int nspin = wf.nspin();
  vector<vector<double> > sp(nspin);

  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    SlaterDet& sd = *(wf.sd(ispin,0));
    const int n = sd.nst();
    if ( erange == 0.0 )
      erange = eVolt * ( sd.eig(n-1) - sd.eig(0) ) + 3 * width;
    const int np = erange / de;
    sp[ispin].resize(np);
    fill(sp[ispin].begin(),sp[ispin].end(),0.0);

    MLWFTransform* mlwft = new MLWFTransform(sd);

    mlwft->update();

    if ( ui->onpe0() )
      cout << "<dipole_matrix_elements ispin=\"" << ispin+1 << "\">" << endl;

    const DoubleMatrix *amat[6];
    for ( int i = 0; i < 6; i++ )
      amat[i] = mlwft->a(i);

    const double *c0 = amat[0]->cvalptr();
    const double *s0 = amat[1]->cvalptr();
    const double *c1 = amat[2]->cvalptr();
    const double *s1 = amat[3]->cvalptr();
    const double *c2 = amat[4]->cvalptr();
    const double *s2 = amat[5]->cvalptr();

    const Context& ctxt = amat[0]->context();

    // loop over global indices (i,j)
    for ( int i = 0; i < n; i++ )
    {
      for ( int j = i+1; j < n; j++ )
      {
        // if i occupied and j empty
        if ( sd.occ(i) > 0.0 && sd.occ(j) == 0.0 )
        {
          const double delta_e = eVolt * ( sd.eig(j) - sd.eig(i) );
          // dipole transition strength w
          double w = 0.0;

          // if element (i,j) is located on the current task,
          // compute local indices (iloc,jloc)
          const int pr = amat[0]->pr(i);
          const int pc = amat[0]->pc(j);
          if ( pr == ctxt.myrow() && pc == ctxt.mycol() )
          {
            const int iloc = amat[0]->l(i) * amat[0]->mb() + amat[0]->x(i);
            const int jloc = amat[0]->m(j) * amat[0]->nb() + amat[0]->y(j);
            const int mloc = amat[0]->mloc();
            // position in local array is k = iloc+mloc*jloc
            const int k = iloc + mloc * jloc;
            const double fac = 0.5 * M_1_PI;
            double c[3] = { fac*c0[k], fac*c1[k], fac*c2[k] };
            double s[3] = { fac*s0[k], fac*s1[k], fac*s2[k] };

            // cc, ss: matrix elements in cartesian coordinates
            double cc[3], ss[3];
            cell.vecmult3x3(cell.amat(),c,cc);
            cell.vecmult3x3(cell.amat(),s,ss);

            w = cc[0]*cc[0]+cc[1]*cc[1]+cc[2]*cc[2]+
                ss[0]*ss[0]+ss[1]*ss[1]+ss[2]*ss[2];

            // add contribution to the absorption spectrum
            for ( int ie = 0; ie < np; ie++ )
            {
              const double t = ( emin + ie * de - delta_e ) / width;
              sp[ispin][ie] += w * width * sqrt(M_PI) * exp(-t*t);
            }

            // only send if not on pe 0
            if ( !ui->onpe0() )
              ctxt.dsend(1,1,&w,1,0,0);

          } // if pr,pc

          // receive if on pe 0 and element was sent from another pe
          if ( ui->onpe0() && !( pr==0 && pc==0 ) )
            ctxt.drecv(1,1,&w,1,pr,pc);

          if ( ui->onpe0() )
          {
            cout.setf(ios::fixed, ios::floatfield);
            cout.setf(ios::right, ios::adjustfield);
            cout << " <dipole i=\"" << i+1 << "\" j=\"" << j+1 << "\"> ";
            cout << setprecision(6) << delta_e << " "
                 << w << " </dipole>" << endl;
          }
        } // if i occupied and j empty
      } // for j
    } // for i

    // sum contributions to sp from all tasks
    ctxt.dsum(np,1,&sp[ispin][0],np);

    if ( ui->onpe0() )
      cout << "</dipole_matrix_elements>" << endl;

    delete mlwft;

  } // ispin

  // write spectrum to file spectrum.dat
  if ( ui->onpe0() )
  {
    ofstream spfile(spfilename);
    for ( int ispin = 0; ispin < nspin; ispin++ )
    {
      spfile << "# spectrum ispin=" << ispin+1
             << " width=" << width << endl;
      for ( int ie = 0; ie < sp[ispin].size(); ie++ )
        spfile << emin + ie * de << " " << sp[ispin][ie] << endl;
      spfile << endl << endl;
    }
    spfile.close();
  }

  return 0;
}
