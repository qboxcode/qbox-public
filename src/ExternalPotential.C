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
// ExternalPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <functional>
using namespace std;

#include "Basis.h"
#include "ExternalPotential.h"
#include "FourierTransform.h"

void ExternalPotential::update(const ChargeDensity& cd)
{
  const Context &cd_ctxt = cd.context();
  MPI_Comm vcomm = cd.vcomm();
  int mype, vcomm_size;
  MPI_Comm_size(vcomm,&vcomm_size);
  MPI_Comm_rank(vcomm,&mype);
  FourierTransform *vft = cd.vft();
  vector<double> vext_read;

  // All tasks with myrow==0 read the potential
  if ( mype == 0 )
  {
    ifstream vfile(filename_.c_str());
    if ( vfile )
    {
      string tmp;
      for ( int i = 0; i < 2; i++ )
        getline(vfile,tmp);  // skip comments

      int natom;
      vfile >> natom;
      getline(vfile,tmp);

      for ( int i = 0; i < 3; i++ )
      {
        vfile >> n_[i];
        getline(vfile,tmp);
      }

      for ( int i = 0; i < natom; i++)
        getline(vfile,tmp);  // skip atom coordinates

      const int n012 = n_[0] * n_[1] * n_[2];
      vext_read.resize(n012);

      // build a min heap to store the largest values of vext
      // to compute its magnitude
      const int heapsize = n012/1000 + 1;  // +1 to avoid getting 0
      vector<double> heap(heapsize, 0.0);

      double tmp_value;
      int counter = 0;
      for ( int nx = 0; nx < n_[0]; nx++ )
        for ( int ny = 0; ny < n_[1]; ny++ )
          for ( int nz = 0; nz < n_[2]; nz++ )
          {
            vfile >> tmp_value;
            if ( counter < heapsize )
              heap[counter] = abs(tmp_value);
            else
            {
              if ( counter == heapsize )
                make_heap(heap.begin(), heap.end());
              if ( abs(tmp_value) > heap[0] )
              {
                pop_heap(heap.begin(), heap.end(), greater<double>());
                heap[heapsize-1] = abs(tmp_value);
                push_heap(heap.begin(), heap.end(), greater<double>());
              }
            }
            const int ir = nx + ny * n_[0] + nz * n_[0] * n_[1];
            vext_read[ir] = tmp_value;
            counter ++;
          }
      vfile.close();
      magnitude_ = accumulate(heap.begin(), heap.end(), 0.0) / heapsize;
    }
    else
    {
      if ( cd_ctxt.onpe0() )
        cout << "  ExternalPotential::update: file not found: "
             << filename_ << endl;
      cd.context().abort(1);
    }
  }
  if ( cd_ctxt.onpe0() )
  {
    cout << "  ExternalPotential::update: read external "
         << "potential from file: " << filename_ << endl;
    cout << "  ExternalPotential::update: grid size "
         << n_[0] << " " << n_[1] << " " << n_[2] << endl;
  }

  // broadcast sizes and magnitude to lower rows
  MPI_Bcast(&n_[0],3,MPI_INT,0,vcomm);
  MPI_Bcast(&magnitude_,1,MPI_DOUBLE,0,vcomm);

  // create a Basis compatible with the vext grid read from file
  const Basis* vbasis = cd.vbasis();
  const UnitCell& cell = vbasis->cell();
  const D3vector b0 = cell.b(0);
  const D3vector b1 = cell.b(1);
  const D3vector b2 = cell.b(2);

  // determine largest ecut compatible with grid and current unit cell
  int n0max = (n_[0]-2)/2;
  int n1max = (n_[1]-2)/2;
  int n2max = (n_[2]-2)/2;
  double ecut0 = 0.5 * norm2(b0) * n0max*n0max;
  double ecut1 = 0.5 * norm2(b1) * n1max*n1max;
  double ecut2 = 0.5 * norm2(b2) * n2max*n2max;
  ecut_ = min(min(ecut0,ecut1),ecut2);

  if ( cd_ctxt.onpe0() )
  {
    cout << "  ExternalPotential::update: magnitude: " << magnitude_ << endl;
    cout << "  ExternalPotential::update: ecut0: " << ecut0 << endl;
    cout << "  ExternalPotential::update: ecut1: " << ecut1 << endl;
    cout << "  ExternalPotential::update: ecut2: " << ecut2 << endl;
    cout << "  ExternalPotential::update: ecut:  " << ecut_ << endl;
  }

  Basis basis(vcomm,D3vector(0,0,0));
  basis.resize(cell,cell,ecut_);
  if ( cd_ctxt.onpe0() )
  {
    cout << "  ExternalPotential::update: np0: " << basis.np(0) << endl;
    cout << "  ExternalPotential::update: np1: " << basis.np(1) << endl;
    cout << "  ExternalPotential::update: np2: " << basis.np(2) << endl;
  }

  assert(basis.np(0)<=vft->np0());
  assert(basis.np(1)<=vft->np1());
  assert(basis.np(2)<=vft->np2());

  // FourierTransform for vext grid
  FourierTransform ft(basis,n_[0],n_[1],n_[2]);

  // scatter parts of vext_read to lower rows
  vector<double> vext_read_loc(ft.np012loc());
  vector<int> scounts(vcomm_size,0);
  vector<int> sdispls(vcomm_size,0);
  int displ = 0;
  for ( int iproc = 0; iproc < vcomm_size; iproc++ )
  {
    sdispls[iproc] = displ;
    scounts[iproc] = ft.np012loc(iproc);
    displ += ft.np012loc(iproc);
  }

  MPI_Scatterv(&vext_read[0],&scounts[0],&sdispls[0],MPI_DOUBLE,
               &vext_read_loc[0],ft.np012loc(),MPI_DOUBLE,0,vcomm);

  vector<complex<double> > tmp_r(ft.np012loc());
  for ( int ir = 0; ir < tmp_r.size(); ir++ )
    tmp_r[ir] = complex<double>(vext_read_loc[ir],0);

  // compute Fourier coefficients
  vector<complex<double> > vext_g_(basis.localsize());
  ft.forward(&tmp_r[0],&vext_g_[0]);
  // vext_g_ now contains the Fourier coefficients of vext

  // interpolate to vft grid
  FourierTransform ft2(basis,vft->np0(),vft->np1(),vft->np2());
  tmp_r.resize(vft->np012loc());
  vext_r_.resize(tmp_r.size());
  ft2.backward(&vext_g_[0],&tmp_r[0]);
  for ( int i = 0; i < vext_r_.size(); i++ )
    vext_r_[i] = real(tmp_r[i]);

  if ( amplitude_ == 0.0 )
  {
    // If amplitude_ = 0.0, use following scheme to get an amplitude.
    // Empirically, an absolute magnitude of 1.0E-3 ~ 1.0E-5 hartree for Vext
    // would be suitable. Here the amplitude is set to scale the
    // magnitude of vext to be 1.0E-3 hartree
    amplitude_ = 1.0E-3 / magnitude_;
    if ( cd_ctxt.onpe0() )
      cout << "  ExternalPotential::update: amplitude automatically determined to be "
           << amplitude_ << endl;
  }
}

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
