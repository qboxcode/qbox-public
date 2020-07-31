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
// WavefunctionHandler.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "WavefunctionHandler.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "Timer.h"
#include "SampleReader.h"
#include "MPIdata.h"

#include "StrX.h"
// XML transcoding for loading grid_functions
#include <xercesc/util/Base64.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;
#include <iostream>
#include <cassert>
#include <cstring> // memcpy
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
WavefunctionHandler::WavefunctionHandler(Wavefunction& wf,
  DoubleMatrix& gfdata, int& current_gfdata_pos ) : wf_(wf), gfdata_(gfdata),
  current_gfdata_pos_(current_gfdata_pos) {}

////////////////////////////////////////////////////////////////////////////////
WavefunctionHandler::~WavefunctionHandler() {}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::byteswap_double(size_t n, double* x)
{
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 )
  {
    unsigned char tmp;
    tmp = c[7]; c[7] = c[0]; c[0] = tmp;
    tmp = c[6]; c[6] = c[1]; c[1] = tmp;
    tmp = c[5]; c[5] = c[2]; c[2] = tmp;
    tmp = c[4]; c[4] = c[3]; c[3] = tmp;

    c+=8;
  }
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::startElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname,
  const Attributes& attributes)
{
  const bool onpe0 = MPIdata::onpe0();
  // cout << " WavefunctionHandler::startElement " << StrX(qname) << endl;
  string locname = StrX(localname).localForm();

  int nspin=1, nel=0, nempty=0;

  // consider only elements that are dealt with directly by WavefunctionHandler

  if ( locname == "wavefunction" || locname == "wavefunction_velocity" )
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      if ( attrname == "ecut")
      {
        ecut = atof(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nspin")
      {
        nspin = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nempty" )
      {
        nempty = atoi(StrX(attributes.getValue(index)).localForm());
      }
      else if ( attrname == "nel" )
      {
        nel = atoi(StrX(attributes.getValue(index)).localForm());
      }
    }

    if ( onpe0 )
      cout << " WavefunctionHandler::startElement: " << locname
           << " nspin=" << nspin << " nel=" << nel << " nempty=" << nempty
           << endl;

    current_ispin = 0;
    current_ikp = 0;

    wf_.set_nel(nel);
    wf_.set_nspin(nspin);
    wf_.set_nempty(nempty);

    // remove default kpoint k=0
    wf_.del_kpoint(D3vector(0.0,0.0,0.0));
  }
  else if ( locname == "domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "a")
      {
        stst >> a;
      }
      else if ( attrname == "b" )
      {
        stst >> b;
      }
      else if ( attrname == "c" )
      {
        stst >> c;
      }
    }

    //cout << " WavefunctionHandler::startElement: domain" << endl;
    uc.set(a,b,c);
    //cout << uc;
  }
  else if ( locname == "reference_domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "a")
      {
        stst >> a;
      }
      else if ( attrname == "b" )
      {
        stst >> b;
      }
      else if ( attrname == "c" )
      {
        stst >> c;
      }
    }

    //cout << " WavefunctionHandler::startElement: reference_domain" << endl;
    ruc.set(a,b,c);
    //cout << ruc;
  }
  else if ( locname == "density_matrix")
  {
    unsigned int len = attributes.getLength();
    int dmat_size = 0;
    string dmat_form;
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "form")
      {
        stst >> dmat_form;
      }
      else if ( attrname == "size" )
      {
        stst >> dmat_size;
      }
    }
    if ( dmat_form != "diagonal" )
    {
      cout << "WavefunctionHandler: density_matrix must be diagonal" << endl;
      MPI_Abort(MPIdata::comm(),1);
    }
    dmat_.resize(dmat_size);
  }
  else if ( locname == "grid")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "nx")
      {
        stst >> nx_;
      }
      else if ( attrname == "ny" )
      {
        stst >> ny_;
      }
      else if ( attrname == "nz" )
      {
        stst >> nz_;
      }
    }

    if ( ecut == 0.0 )
    {
      // ecut attribute was not specified. Infer from grid size
      // Ecut = max(ecut_x,ecut_y,ecut_z)

      // When importing grids with Dirichlet b.c. grid sizes can be odd
      // round nx,ny,nz to next even number to compute ecut
      // use nx+nx%2 instead of nx
      double g0_max = ((2*(nx_+nx_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(0));
      double g1_max = ((2*(ny_+ny_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(1));
      double g2_max = ((2*(nz_+nz_%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(2));
      double ecut0 = 0.125 * g0_max * g0_max;
      double ecut1 = 0.125 * g1_max * g1_max;
      double ecut2 = 0.125 * g2_max * g2_max;

      ecut = max(max(ecut0,ecut1),ecut2);
      cout << " ecut=" << 2*ecut << " Ry" << endl;
    }

    wf_.resize(uc,ruc,ecut);
  }
  else if ( locname == "slater_determinant")
  {
    if ( onpe0 )
      cout << " WavefunctionHandler::startElement: slater_determinant" << endl;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "kpoint")
      {
        stst >> current_kx >> current_ky >> current_kz;
        if ( onpe0 )
          cout << " kpoint=" << current_kx
               << " " << current_ky << " " << current_kz;
      }
      else if ( attrname == "weight" )
      {
        stst >> current_weight;
        if ( onpe0 )
          cout << " weight=" << current_weight;
      }
      else if ( attrname == "size" )
      {
        last_size = current_size;
        stst >> current_size;
        if ( onpe0 )
          cout << " size=" << current_size;
      }
      else if ( attrname == "spin" )
      {
        std::string spin;
        stst >> spin;
        if (spin == "up" )
        {
          current_ispin=0;
        }
        else if (spin == "down" )
        {
          assert(wf_.nspin() == 2);
          int last_ispin=current_ispin;
          current_ispin=1;
          // reset kpoint index if necessary
          if ( last_ispin==0 && current_ispin==1 )
          {
            current_ikp=0;
          }
        }
        if ( onpe0 )
          cout << " read slater_determinant: spin=" << spin;
      }
    }

    if ( onpe0 )
      cout << endl;

    // All SlaterDets of a given spin must have the same size
    // check if current_size differs from previous size
    if ( ( current_ikp != 0 ) && ( current_size != last_size ) )
    {
      cout << "SlaterDet size differs from previous size" << endl;
      MPI_Abort(MPIdata::comm(),1);
    }

    // add kpoint only if spin is up (if spin is down,
    // the kpoint should be defined already)
    if ( current_ispin == 0 )
      wf_.add_kpoint(D3vector(current_kx,current_ky,current_kz),current_weight);

    // resize SlaterDet if local
    int isp_loc = wf_.isp_local(current_ispin);
    int ikp_loc = wf_.ikp_local(current_ikp);
    if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
    {
      assert(wf_.sd(isp_loc,ikp_loc));
      wf_.sd(isp_loc,ikp_loc)->resize(wf_.cell(),
             wf_.refcell(), wf_.ecut(), current_size);
    }

    wf_.nst_[current_ispin] = current_size;
  }
  else if ( locname == "grid_function")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname = StrX(attributes.getLocalName(index)).localForm();
      string attrval = StrX(attributes.getValue(index)).localForm();
      istringstream stst(attrval);
      if ( attrname == "nx")
      {
        stst >> current_gf_nx;
      }
      else if ( attrname == "ny" )
      {
        stst >> current_gf_ny;
      }
      else if ( attrname == "nz" )
      {
        stst >> current_gf_nz;
      }
      else if ( attrname == "encoding" )
      {
        stst >> current_gf_encoding;
      }
    }

    //cout << " WavefunctionHandler::startElement: grid_function"
    //     << " nx=" << current_gf_nx
    //     << " ny=" << current_gf_ny
    //     << " nz=" << current_gf_nz
    //     << "\n encoding=" << current_gf_encoding
    //     << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::endElement(const XMLCh* const uri,
  const XMLCh* const localname, const XMLCh* const qname, string& content)
{
  const Context& sd_ctxt = wf_.sd_context();
  bool onpe0 = MPIdata::onpe0();
  string locname = StrX(localname).localForm();
  //cout << " WavefunctionHandler::endElement " << locname << endl;
  if ( locname == "density_matrix")
  {
    istringstream stst(content);
    for ( int i = 0; i < dmat_.size(); i++ )
      stst >> dmat_[i];
  }
  else if ( locname == "grid_function")
  {
    // current implementation accepts only full grids as declared in
    // the wavefunction <grid> element
    assert(current_gf_nx==nx_ &&
           current_gf_ny==ny_ &&
           current_gf_nz==nz_ );
  }
  else if ( locname == "slater_determinant")
  {
    int isp_loc = wf_.isp_local(current_ispin);
    int ikp_loc = wf_.ikp_local(current_ikp);
    if ( ( isp_loc >= 0 ) && ( ikp_loc >= 0 ) )
      wf_.sd(isp_loc,ikp_loc)->set_occ(dmat_);
    if ( onpe0 )
      cout << " WavefunctionHandler::endElement: slater_determinant" << endl;
    current_ikp++;
  }
  else if ( locname == "wavefunction" || locname == "wavefunction_velocity" )
  {
    // copy data from gfdata_ to locat wftmpr
    vector<double> wftmpr(gfdata_.mloc());
    // jsrc: column of gfdata_
    int jsrc = 0;
    for ( int ispin = 0; ispin < wf_.nspin(); ++ispin )
    {
      const int nst = wf_.nst(ispin);
      // nb: size of blocks in sd
      const int nb = nst/MPIdata::nstb() + (nst%MPIdata::nstb() > 0 ? 1 : 0);
      int isp_loc = wf_.isp_local(ispin);
      for ( int ikp = 0; ikp < wf_.nkp(); ++ikp )
      {
        MPI_Barrier(MPI_COMM_WORLD);
        int ikp_loc = wf_.ikp_local(ikp);

        SlaterDet* sd = 0;
        FourierTransform* ft = 0;
        vector<complex<double> > wftmp;
        if ( ( isp_loc >= 0 ) && ( ikp_loc ) >= 0 )
        {
          assert(wf_.sd(isp_loc,ikp_loc));
          sd = wf_.sd(isp_loc,ikp_loc);
          //cout << " sd->c().size() " << sd->c().size() << endl;
          const Basis& basis = sd->basis();
          ft = new FourierTransform(basis,nx_,ny_,nz_);
          wftmp.resize((ft->np012loc()));
          //cout << MPIdata::rank() << ": ft->np012loc()="
          //     << ft->np012loc() << endl;
        }

        for ( int n = 0; n < wf_.nst(ispin); ++n )
        {
          MPI_Barrier(MPI_COMM_WORLD);

          // compute MPI rank of sender
          // process column holding column jsrc on the gfdata context
          int csrc = gfdata_.pc(jsrc);
          // n_loc: local column index of jsrc on gfdata
          int n_loc = gfdata_.y(jsrc);
          // src_rank: MPI rank in MPIdata::comm()
          int src_rank = MPIdata::igb() + csrc * MPIdata::ngb();
          assert(src_rank >= 0);
          assert(src_rank < MPIdata::size());

          // compute MPI rank of receiver
          // sd_col: process col of the proc holding state n in the sd context
          int sd_col = n / nb;
          int dest_rank = MPIdata::igb() +
                          MPIdata::ngb()  * ( sd_col +
                          MPIdata::nstb() * ( ikp % MPIdata::nkpb() +
                          MPIdata::nkpb() *   ispin % MPIdata::nspb() ) );
          assert(dest_rank >= 0);
          assert(dest_rank < MPIdata::size());

          int len = gfdata_.mloc();
          bool iamsending = ( src_rank == MPIdata::rank() );
          bool iamreceiving = ( dest_rank == MPIdata::rank() );
          if ( iamsending && iamreceiving )
          {
#ifdef DEBUG
            cout << MPIdata::rank() << ": n=" << n
                 << " jsrc=" << jsrc << " local copy" << endl;
#endif
            memcpy(&wftmpr[0],gfdata_.valptr(n_loc*len),len*sizeof(double));
          }
          else
          {
            if ( iamreceiving )
            {
#ifdef DEBUG
              cout << MPIdata::rank() << ": n=" << n
                   << " jsrc=" << jsrc
                   << " receiving from " << src_rank << endl;
#endif
              MPI_Status stat;
              MPI_Recv(&wftmpr[0],len,MPI_DOUBLE,src_rank,src_rank,
                       MPIdata::comm(),&stat);
            }
            if ( iamsending )
            {
#ifdef DEBUG
              cout << MPIdata::rank() << ": n=" << n
                   << " jsrc=" << jsrc
                   << " sending to " << dest_rank << endl;
#endif
              assert(n_loc < gfdata_.nloc());
              MPI_Send(gfdata_.valptr(n_loc*len),len,MPI_DOUBLE,dest_rank,
                 src_rank,MPIdata::comm());
            }
          }

          // transform state in wftmpr to Fourier coefficients
          if ( sd )
          {
            const Basis& basis = sd->basis();

            // check if n is local to this task
            if ( ( n / nb ) == sd_ctxt.mycol() )
            {
              // copy column of wftmpr to complex array wftmp
              if ( basis.real() )
              {
                // function is real and must be expanded
                for ( int i = 0; i < ft->np012loc(); i++ )
                  wftmp[i] = wftmpr[i];
              }
              else
              {
                // function is complex
                double* p = &wftmpr[0];
                for ( int i = 0; i < ft->np012loc(); i++ )
                  wftmp[i] = complex<double>(p[2*i],p[2*i+1]);
              }
              ComplexMatrix& c = sd->c();
              int nloc = c.y(n);
              ft->forward(&wftmp[0],c.valptr(c.mloc()*nloc));
            }
          }
          jsrc++;
        } // n
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
        cout << MPIdata::rank() << ": end copying ispin=" << ispin
             << " ikp=" << ikp << endl;
#endif
        delete ft;
      } // ikp
    } // ispin


    // if nspin==2, adjust deltaspin_ to reflect the number of states
    // of each spin that were read
    if ( wf_.nspin() == 2 )
    {
      // assume that up spin is the majority spin
      assert(wf_.nst(0) >= wf_.nst(1));
      // The following line is correct for both
      // even and odd number of electrons
      wf_.deltaspin_ = (wf_.nst(0)-wf_.nst(1))/2;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
StructureHandler* WavefunctionHandler::startSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const Attributes& attributes)
{
  // check if element qname can be processed by another StructureHandler
  // If it can, return a pointer to the StructureHandler, otherwise return 0
  // cout << " WavefunctionHandler::startSubHandler " << StrX(qname) << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void WavefunctionHandler::endSubHandler(const XMLCh* const uri,
    const XMLCh* const localname, const XMLCh* const qname,
    const StructureHandler* const last)
{
  string locname = StrX(localname).localForm();
  //cout << " WavefunctionHandler::endSubHandler " << locname << endl;
  delete last;
}
