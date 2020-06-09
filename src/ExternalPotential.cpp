////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2016 The Regents of the University of California
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
// ExternalPotential.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

#include "Basis.h"
#include "ExternalPotential.h"
#include "FourierTransform.h"
#include "Function3d.h"
#include "Base64Transcoder.h"

////////////////////////////////////////////////////////////////////////////////
bool abs_compare(const double &a, const double &b)
{
  return (abs(a) < abs(b));
}

////////////////////////////////////////////////////////////////////////////////
void ExternalPotential::update(const ChargeDensity& cd)
{
  const Context* ctxt = s_.wf.spincontext();
  bool onpe0 = ctxt->onpe0();
  int nprow = ctxt->nprow();
  int myrow = ctxt->myrow();
  int mycol = ctxt->mycol();
  MPI_Comm vcomm = cd.vcomm();

  Timer tm_read_vext;
  double time, tmin, tmax;

  // In cube mode and xml mode, the external potential is
  // read by processors on row 0 and stored in vext_read
  vector<double> vext_read, vext_read_loc;

  tm_read_vext.start();
  if ( fmt_ == "cube" )
  {
    // read cube file, n_'s are determined by cube file
    if ( myrow == 0 )
    {
      ifstream vfile(filename_.c_str());
      if (!vfile)
      {
        if (mycol == 0)
          cout << "  ExternalPotential::update: file not found: "
               << filename_ << endl;
        ctxt->abort(1);
      }
      string tmpstr;
      for (int i = 0; i < 2; i++)
        getline(vfile, tmpstr);  // skip comments
      int natom;
      vfile >> natom;
      getline(vfile, tmpstr);
      for (int i = 0; i < 3; i++)
      {
        vfile >> n_[i];
        getline(vfile, tmpstr);
      }
      for (int i = 0; i < natom; i++)
        getline(vfile, tmpstr);  // skip atom coordinates
      vext_read.resize(n_[0] * n_[1] * n_[2]);
      for (int nx = 0; nx < n_[0]; nx++)
        for (int ny = 0; ny < n_[1]; ny++)
          for (int nz = 0; nz < n_[2]; nz++)
          {
            const int ir = nx + ny * n_[0] + nz * n_[0] * n_[1];
            vfile >> vext_read[ir];
          }
      vfile.close();
    }
    MPI_Bcast(&n_[0],3,MPI_INT,0,vcomm);
  }
  else if ( fmt_ == "xml" )
  {
    if (myrow == 0)
    {
      Function3d f;
      f.read(filename_);
      vext_read = f.val;
      n_[0] = f.nx;
      n_[1] = f.ny;
      n_[2] = f.nz;
    }
    MPI_Bcast(&n_[0],3,MPI_INT,0,vcomm);
  }
  tm_read_vext.stop();
  // at this point, all processes have correct n_ regardless of io mode

  // now construct 2 Fourier transforms ft1 and ft2
  // ft1 is used to transform vext_read_loc to vext_g (G space)
  // ft2 is used to transform vext_g to vext_r_ (R space)
  // the whole process is a Fourier interpolation/extrapolation
  // of the external potential on the charge density grid

  // create a Basis with largest possible ecut that is compatible with
  // the external potential grid from file
  const UnitCell& cell = cd.vbasis()->cell();
  int n0max = (n_[0]-2)/2;
  int n1max = (n_[1]-2)/2;
  int n2max = (n_[2]-2)/2;
  double ecut0 = 0.5 * norm2(cell.b(0)) * n0max*n0max;
  double ecut1 = 0.5 * norm2(cell.b(1)) * n1max*n1max;
  double ecut2 = 0.5 * norm2(cell.b(2)) * n2max*n2max;
  ecut_ = min(min(ecut0,ecut1),ecut2);
  Basis basis(vcomm,D3vector(0,0,0));
  basis.resize(cell,cell,ecut_);

  FourierTransform *vft = cd.vft();

  FourierTransform ft1(basis,n_[0],n_[1],n_[2]);
  vext_read_loc.resize(ft1.np012loc());
  vector<complex<double> > vext_g(basis.localsize());

  // check that the basis fits in the vft grid
  //assert(basis.fits_in_grid(vft->np0(),vft->np1(),vft->np2()));

  FourierTransform ft2(basis,vft->np0(),vft->np1(),vft->np2());
  vext_r_.resize(ft2.np012loc());

  // xml mode or cube mode: processors on row 0 scatter
  // vext to other rows
  vector<int> scounts(nprow,0);
  vector<int> sdispls(nprow,0);
  int displ = 0;
  for ( int iproc = 0; iproc < nprow; iproc++ )
  {
    sdispls[iproc] = displ;
    scounts[iproc] = ft1.np012loc(iproc);
    displ += ft1.np012loc(iproc);
  }
  MPI_Scatterv(&vext_read[0],&scounts[0],&sdispls[0],MPI_DOUBLE,
               &vext_read_loc[0],ft1.np012loc(),MPI_DOUBLE,0,vcomm);

  // now vext_read_loc on all processors contains the correct portion of vext
  // Fourier forward transform vext_read_loc to vext_g
  vector<complex<double> > tmp_r(ft1.np012loc());
  for ( int ir = 0; ir < tmp_r.size(); ir++ )
    tmp_r[ir] = complex<double>(vext_read_loc[ir],0);
  ft1.forward(&tmp_r[0],&vext_g[0]);
  // Fourier backward transform vext_g to vext_r_
  tmp_r.resize(ft2.np012loc());
  ft2.backward(&vext_g[0],&tmp_r[0]);
  for ( int i = 0; i < vext_r_.size(); i++ )
    vext_r_[i] = real(tmp_r[i]);
  if ( onpe0 )
  {
    cout << "  ExternalPotential::update: read external potential from file:  "
         << filename_ << endl;
    cout << "  ExternalPotential::update: grid size n0 = "
         << n_[0] << ", n1 = " << n_[1] << ", n2 = " << n_[2] << endl;
    cout << "  ExternalPotential::update: ecut:  " << ecut_ << endl;
    if ( amplitude_ != 0 )
      cout << "  ExternalPotential::update: amplitude:  " << amplitude_ << endl;
  }
  if ( amplitude_ == 0.0 )
  {
    // If amplitude_ = 0.0, use following scheme to get an amplitude.
    // Empirically, an absolute magnitude of 1.0E-3 ~ 1.0E-5 hartree for Vext
    // would be suitable. Here the amplitude is set to scale the
    // magnitude of vext to be 1.0E-3 hartree
    if (vext_read_loc.size() > 0)
      magnitude_ = abs(*max_element(vext_read_loc.begin(),
                       vext_read_loc.end(), abs_compare));
    ctxt->dmax('C',1,1,&magnitude_,1);
    MPI_Bcast(&magnitude_,1,MPI_DOUBLE,0,vcomm);
    amplitude_ = 1.0E-3 / magnitude_;
    if ( onpe0 )
      cout << "  ExternalPotential::update: amplitude automatically"
           << " determined to be " << amplitude_
           << " (max abs value of vext = " << magnitude_ << ")" << endl;
  }

  time = tm_read_vext.real();
  tmin = time;
  tmax = time;
  ctxt->dmin(1,1,&tmin,1);
  ctxt->dmax(1,1,&tmax,1);
  if ( onpe0 )
  {
    cout << "  ExternalPotential::update: vext read time "
         << "min: " << tmin << " max: " << tmax << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
double ExternalPotential::compute_eext(const ChargeDensity& cd)
{
  // Eext =  integral ( rhor * vext_r )
  double eext = 0.0;
  double omega = s_.wf.cell().volume();
  const int np012loc = cd.vft()->np012loc();
  const int np012 = cd.vft()->np012();
  const std::vector<std::vector<double> > &rhor = cd.rhor;
  for ( int ispin = 0; ispin < s_.wf.nspin(); ispin++ )
  {
    for ( int ir = 0; ir < np012loc; ir++ )
    {
      assert( rhor[ispin].size() == np012loc );
      assert( vext_r_.size() == np012loc );
      eext += rhor[ispin][ir] * v(ir);
    }
  }
  double eext_sum = 0.0;
  MPI_Allreduce(&eext,&eext_sum,1,MPI_DOUBLE,MPI_SUM,cd.vbasis()->comm());
  eext_sum =  eext_sum * omega / np012;
  return eext_sum;
}
