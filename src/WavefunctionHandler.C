////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionHandler.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionHandler.C,v 1.10 2004-09-14 22:24:11 fgygi Exp $

#if USE_XERCES

#include "WavefunctionHandler.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "Timer.h"

#include "StrX.h"
// XML transcoding for loading grid_functions
#include <xercesc/util/Base64.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;
#include <iostream>
#include <cassert>
#include <sstream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
WavefunctionHandler::WavefunctionHandler(Wavefunction& wf, 
  DoubleMatrix& gfdata) : wf_(wf), gfdata_(gfdata), ecut(0.0)
{
  // if the gfdata matrix has finite dimensions, set read_from_gfdata flag
  read_from_gfdata = ( gfdata.n() > 0 );
  current_igf = 0;
}

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
  // cout << " WavefunctionHandler::startElement " << StrX(qname) << endl;
  string locname(XMLString::transcode(localname));
  
  int nspin=1, nel=0, nempty=0;
  // consider only elements that are dealt with directly by WavefunctionHandler
  
  if ( locname == "wavefunction" || locname == "wavefunction_velocity" )
  {
    
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
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
    
    cout << " WavefunctionHandler::startElement: " << locname
         << " nspin=" << nspin << " nel=" << nel << " nempty=" << nempty
         << endl;
    
    current_ispin = 0;
    current_ikp = 0;
    current_n = 0;
    
    wf_.set_nel(nel);
    wf_.set_nspin(nspin);
    wf_.set_nempty(nempty);
    
    int tag;
    // notify listening nodes
    if ( locname == "wavefunction" )
      tag = 3;
    else
      tag = 4; // wavefunction_velocity
    
    wf_.context().ibcast_send(1,1,&tag,1);
    wf_.context().ibcast_send(1,1,&nel,1);
    wf_.context().ibcast_send(1,1,&nspin,1);
    wf_.context().ibcast_send(1,1,&nempty,1);
    
    // current implementation for nspin = 1 and nkpoint = 1
    assert(nspin==1);
  }
  else if ( locname == "domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
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
    
    // notify listening nodes
    double buf[9];
    buf[0] = uc.a(0).x; buf[1] = uc.a(0).y; buf[2] = uc.a(0).z;
    buf[3] = uc.a(1).x; buf[4] = uc.a(1).y; buf[5] = uc.a(1).z;
    buf[6] = uc.a(2).x; buf[7] = uc.a(2).y; buf[8] = uc.a(2).z;
    wf_.context().dbcast_send(9,1,buf,1);
  }
  else if ( locname == "reference_domain")
  {
    D3vector a,b,c;
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
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
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
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
      wf_.context().abort(1);
    }
    dmat.resize(dmat_size);
  }
  else if ( locname == "grid")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
      istringstream stst(attrval);
      if ( attrname == "nx")
      {
        stst >> nx;
      }
      else if ( attrname == "ny" )
      {
        stst >> ny;
      }
      else if ( attrname == "nz" )
      {
        stst >> nz;
      }
    }
    
    if ( ecut == 0.0 )
    {
      // ecut attribute was not specified. Infer from grid size
      // Ecut = max(ecut_x,ecut_y,ecut_z)
    
      // When importing grids with Dirichlet b.c. grid sizes can be odd
      // round nx,ny,nz to next even number to compute ecut
      // use nx+nx%2 instead of nx
      double g0_max = ((2*(nx+nx%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(0));
      double g1_max = ((2*(ny+ny%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(1));
      double g2_max = ((2*(nz+nz%2)-1.0)/2.0) * 2.0 * M_PI / length(uc.a(2));
      double ecut0 = 0.125 * g0_max * g0_max;
      double ecut1 = 0.125 * g1_max * g1_max;
      double ecut2 = 0.125 * g2_max * g2_max;

      ecut = max(max(ecut0,ecut1),ecut2);
      cout << " ecut=" << 2*ecut << " Ry" << endl;
    }
    
    // notify listening nodes of ecut
    wf_.context().dbcast_send(1,1,&ecut,1);
    
    // notify listening nodes of the reference_domain
    // note: the reference_domain is optional in the sample file
    // notify listening nodes
    double buf[9];
    buf[0] = ruc.a(0).x; buf[1] = ruc.a(0).y; buf[2] = ruc.a(0).z;
    buf[3] = ruc.a(1).x; buf[4] = ruc.a(1).y; buf[5] = ruc.a(1).z;
    buf[6] = ruc.a(2).x; buf[7] = ruc.a(2).y; buf[8] = ruc.a(2).z;
    wf_.context().dbcast_send(9,1,buf,1);
    
    wf_.resize(uc,ruc,ecut);
  }
  else if ( locname == "slater_determinant")
  {
    const Basis& basis = wf_.sd(current_ispin,current_ikp)->basis();
    // check the size of the basis generated
    //cout << " sd basis: np0,np1,np2 = " << basis.np(0)
    //     << " " << basis.np(1)
    //     << " " << basis.np(2)
    //     << endl;
    ft = new FourierTransform(basis,basis.np(0),basis.np(1),basis.np(2));
    wftmp.resize((ft->np012loc()));
  }
  else if ( locname == "grid_function")
  {
    unsigned int len = attributes.getLength();
    for (unsigned int index = 0; index < len; index++)
    {
      string attrname(XMLString::transcode(attributes.getLocalName(index)));
      string attrval(XMLString::transcode(attributes.getValue(index)));
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
  string locname(XMLString::transcode(localname));
  //cout << " WavefunctionHandler::endElement " << locname << endl;
  if ( locname == "density_matrix")
  {
    istringstream stst(content);
    for ( int i = 0; i < dmat.size(); i++ )
      stst >> dmat[i];
      
    // send dmat to listening nodes
    //!! this works only for 1 kpoint, 1 spin
    SlaterDet* sd = wf_.sd(current_ispin,current_ikp);
    assert(sd != 0);
    sd->context().dbcast_send(dmat.size(),1,&dmat[0],1);
    sd->set_occ(dmat);
  }
  else if ( locname == "grid_function")
  {
    // current implementation accepts only full grids
    assert(current_gf_nx==ft->np0());
    assert(current_gf_ny==ft->np1());
    assert(current_gf_nz==ft->np2());
    
    if ( read_from_gfdata )
    {
      // do nothing
      //cout << "WavefunctionHandler::endElement: current_igf=" << current_igf
      //     << endl;
      current_igf++;
    }
    else
    {
      if ( current_gf_encoding == "text" )
      {
        istringstream stst(content);
        valarray<double> wftmpr(current_gf_nx*current_gf_ny*current_gf_nz);
        int ii = 0;
        for ( int k = 0; k < current_gf_nz; k++ )
          for ( int j = 0; j < current_gf_ny; j++ )
            for ( int i = 0; i < current_gf_nx; i++ )
            {
              stst >> wftmpr[ii++];
            }
        // send subgrids to listening nodes
 
        SlaterDet* sd = wf_.sd(current_ispin,current_ikp);
        assert(sd != 0);
        // pcol = process column destination
        int pcol = sd->c().pc(current_n);
        for ( int prow = 0; prow < sd->context().nprow(); prow++ )
        {
          int size = ft->np2_loc(prow) * ft->np0() * ft->np1();
          int istart = ft->np2_first(prow) * ft->np0() * ft->np1();
          // send subgrid to node (prow,pcol)
          if ( !(prow==0 && pcol==0) )
          {
            //cout << sd->context();
            //cout << sd->context().mype() << ": sending subgrid size to process "
            //     << "(" << prow << "," << pcol << ")" << endl;
            sd->context().isend(1,1,&size,1,prow,pcol);
 
            //cout << sd->context().mype() << ": sending subgrid to process "
            //     << "(" << prow << "," << pcol << ")" << endl;
            sd->context().dsend(size,1,&wftmpr[istart],1,prow,pcol);
            //cout << sd->context().mype() << ": subgrid sent "
            //     << endl;
          }
        }
 
        // if destination column is pcol=0, copy to complex array on node 0
        // and process
        if ( pcol == 0 )
        {
          for ( int i = 0; i < ft->np012loc(); i++ )
          {
            wftmp[i] = wftmpr[i];
          }
          ComplexMatrix& c = wf_.sd(current_ispin,current_ikp)->c();
          ft->forward(&wftmp[0],c.valptr(c.mloc()*current_n));
        }
      }
      else if ( current_gf_encoding == "base64" )
      {
        unsigned int length;
        Timer tm;
        tm.start();
        XMLByte* b = Base64::decode((XMLByte*)content.c_str(), &length);
        tm.stop();
        // cout << " Base64::decode time: " << tm.real() << endl;
        assert(b!=0);
        // cout << " base64 segment length: " << length << endl;
 
        // use data in b
        assert(length/sizeof(double)==ft->np012());
 
        double* d = (double*) b;
#if AIX
        tm.reset();
        tm.start();
        byteswap_double(ft->np012(),d);
        tm.stop();
        // cout << " byteswap time: " << tm.real() << endl;
#endif

        SlaterDet* sd = wf_.sd(current_ispin,current_ikp);
        assert(sd != 0);
        // pcol = process column destination
        tm.reset();
        tm.start();
        int pcol = sd->c().pc(current_n);
        for ( int prow = 0; prow < sd->context().nprow(); prow++ )
        {
          int size = ft->np2_loc(prow) * ft->np0() * ft->np1();
          int istart = ft->np2_first(prow) * ft->np0() * ft->np1();
          // send subgrid to node (prow,pcol)
          if ( !(prow==0 && pcol==0) )
          {
            //cout << sd->context();
            //cout << sd->context().mype() << ": sending subgrid size to process "
            //     << "(" << prow << "," << pcol << ")" << endl;
            sd->context().isend(1,1,&size,1,prow,pcol);
 
            //cout << sd->context().mype() << ": sending subgrid to process "
            //     << "(" << prow << "," << pcol << ")" << endl;
            sd->context().dsend(size,1,&d[istart],1,prow,pcol);
            //cout << sd->context().mype() << ": subgrid sent "
            //     << endl;
          }
        }
        tm.stop();
        // cout << " send time: " << tm.real() << endl;
 
        // if destination column is pcol=0, copy to complex array on node 0
        // and process
        if ( pcol == 0 )
        {
          tm.reset();
          tm.start();
          for ( int i = 0; i < ft->np012loc(); i++ )
          {
            wftmp[i] = d[i];
          }
          XMLString::release(&b);
          ComplexMatrix& c = wf_.sd(current_ispin,current_ikp)->c();
          ft->forward(&wftmp[0],c.valptr(c.mloc()*current_n));
          tm.stop();
          // cout << " process time: " << tm.real() << endl;
        }
      }
      else
      {
        cout << "WavefunctionHandler: unknown grid_function encoding" << endl;
        return;
      }
    }
    //cout << " grid_function read n=" << current_n << endl;
    current_n++;
  }
  else if ( locname == "slater_determinant")
  {
    delete ft;
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
  string locname(XMLString::transcode(localname));
  //cout << " WavefunctionHandler::endSubHandler " << locname << endl;
  delete last;
}
      
#endif
