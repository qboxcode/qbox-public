
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
using namespace std;

#include "Basis.h"
#include "ExternalPotential.h"
#include "FourierTransform.h"

void readRealSpaceFunction(const Context& ctxt, vector<double>& fr,
                           const FourierTransform& ft,
                           const string filename, const bool parallel);

void ExternalPotential::update(const ChargeDensity& cd)
{
  const Context& ctxt = *(s_.wf.spincontext());

  // process 0 reads the grid size and bcast to all processes
  if ( ctxt.onpe0() )
  {
    ifstream vfile(filename_.c_str());
    if ( vfile )
      cout << "  ExternalPotential::update: read external "
           << "potential from file: " << filename_ << endl;
    else
    {
      cout << "  ExternalPotential::update: file not found: "
           << filename_ << endl;
      ctxt.abort(1);
    }
    vfile >> n_[0] >> n_[1] >> n_[2];
    cout << "  ExternalPotential::update: grid size "
         << n_[0] << " " << n_[1] << " " << n_[2] << endl;
    vfile.close();
  }
  MPI_Bcast(&n_[0],3,MPI_INT,0,cd.vcomm());

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

  if ( ctxt.onpe0() )
  {
    cout << "  ExternalPotential::update: ecut0: " << ecut0 << endl;
    cout << "  ExternalPotential::update: ecut1: " << ecut1 << endl;
    cout << "  ExternalPotential::update: ecut2: " << ecut2 << endl;
    cout << "  ExternalPotential::update: ecut:  " << ecut_ << endl;
  }

  Basis basis(cd.vcomm(),D3vector(0,0,0));
  basis.resize(cell,cell,ecut_);
  if ( ctxt.onpe0() )
  {
    cout << "  ExternalPotential::update: np0: " << basis.np(0) << endl;
    cout << "  ExternalPotential::update: np1: " << basis.np(1) << endl;
    cout << "  ExternalPotential::update: np2: " << basis.np(2) << endl;
  }

  // build FT to transform vext to g space
  FourierTransform ft(basis,n_[0],n_[1],n_[2]);

  // read vext from file

  vector<double> vext_r_read;
  bool parallel_read = false; // later parallel_read should be determined by vext file
  readRealSpaceFunction(ctxt,vext_r_read,ft,filename_, parallel_read);

  vector<complex<double> > tmp_r(ft.np012loc());
  for ( int ir = 0; ir < tmp_r.size(); ir++ )
    tmp_r[ir] = complex<double>(vext_r_read[ir],0);

  // compute Fourier coefficients
  vector<complex<double> > vext_g_(basis.localsize());
  ft.forward(&tmp_r[0],&vext_g_[0]);
  // vext_g_ now contains the Fourier coefficients of vext

  // interpolate to vft grid of charge density
  FourierTransform *vft = cd.vft();
  FourierTransform ft2(basis,vft->np0(),vft->np1(),vft->np2());
  tmp_r.resize(vft->np012loc());
  vext_r_.resize(tmp_r.size());
  ft2.backward(&vext_g_[0],&tmp_r[0]);
  for ( int i = 0; i < vext_r_.size(); i++ )
    vext_r_[i] = real(tmp_r[i]);
}

void readRealSpaceFunction(const Context& ctxt, vector<double>& fr,
                           const FourierTransform& ft,
                           const string filename, const bool parallel)
{
  // This function read a real space function from a file and distribute
  // it to fr variable among all process of the context
  // ft is used to determine how fr is stored at each process
  // If parallel is true, file will be read in parallel assuming base64 encoding,
  // otherwise file will be read by process 0 in text format,
  // and then distribute to other processes

  // generate column and row communicators
  MPI_Comm comm = ctxt.comm();
  MPI_Comm comm_col;
  MPI_Comm comm_row;
  MPI_Comm_split(comm,ctxt.mycol(),ctxt.myrow(),&comm_col);
  MPI_Comm_split(comm,ctxt.myrow(),ctxt.mycol(),&comm_row);

  if ( ctxt.onpe0())
  {
    assert( ctxt.mycol() == 0 );
    fr.resize(ft.np012());
  }
  else
    fr.resize(ft.np012loc());

  if ( parallel )
  {
    if ( ctxt.mycol() == 0 )
    {
      MPI_Info info;
      MPI_File vfile;
      int err;
      err = MPI_File_open(comm_col,filename.c_str(),MPI_MODE_RDONLY,info,&vfile);
      if ( err != 0 )
      {
        cout << " Error opening file " << filename << endl;
        ctxt.abort(1);
      }

      // todo: implement MPI read

      MPI_File_close(&vfile);
    }
  } // if parallel
  else
  {
    if ( ctxt.onpe0() )
    {
      ifstream vfile(filename.c_str());
      int n1, n2, n3;
      vfile >> n1 >> n2 >> n3;
      for ( int ir = 0; ir < fr.size(); ir++ )
        vfile >> fr[ir];
      vfile.close();
    }

    void* sbuf;
    void* rbuf;
    if ( ctxt.onpe0() )
    {
      sbuf = &fr[0];
      rbuf = MPI_IN_PLACE;
    }
    else
      rbuf = &fr[0];

    int displ = 0;
    vector<int> scounts(ctxt.nprow(),0);
    vector<int> displs(ctxt.nprow(),0);
    for ( int iproc = 0; iproc < ctxt.nprow(); iproc++ )
    {
      displs[iproc] = displ;
      displ += ft.np012loc(iproc);
      scounts[iproc] = ft.np012loc(iproc);
    }

    // process 0 scatter fr to other processes on first column
    MPI_Barrier(comm);
    if ( ctxt.mycol() == 0 )
    {
      MPI_Scatterv(sbuf,&scounts[0],&displs[0],MPI_DOUBLE,
                   rbuf,ft.np012loc(),MPI_DOUBLE,0,comm_col);
    }
    if ( ctxt.onpe0() )
    {
      fr.resize(ft.np012loc());
    }

    // first column bcast fr to other columns
    MPI_Bcast(&fr[0],ft.np012loc(),MPI_DOUBLE,0,comm_row);
  } // if parallel
}
