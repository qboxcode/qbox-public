////////////////////////////////////////////////////////////////////////////////
//
// SampleReader.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleReader.C,v 1.11 2003-12-02 20:26:05 fgygi Exp $


#include "Sample.h"
#include "SampleReader.h"
#include "SpeciesReader.h"
#include "StructuredDocumentHandler.h"
#include "SampleHandler.h"
#include "Basis.h"
#include "FourierTransform.h"
#include "SlaterDet.h"

#include "XMLGFPreprocessor.h"

#include "Timer.h"

#include <cassert>
#include <string>
#include <iostream>
#include <sys/stat.h>
using namespace std;

#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
using namespace xercesc;

////////////////////////////////////////////////////////////////////////////////
SampleReader::SampleReader(const Context& ctxt) : ctxt_(ctxt) {}

////////////////////////////////////////////////////////////////////////////////
void SampleReader::readSample (Sample& s, const string uri, bool serial)
{
  const char* encodingName = "UTF-8";
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Auto;
  SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Always;
  //SAX2XMLReader::ValSchemes valScheme = SAX2XMLReader::Val_Never;
  bool expandNamespaces = false;
  bool doNamespaces = true;
  bool doSchema = true;
  bool schemaFullChecking = true;
  bool namespacePrefixes = false;
  SAX2XMLReader* parser = 0;
  
  Wavefunction wfvtmp(s.wf);
  
  MemBufInputSource* memBufIS = 0;
  
  // XMLPlatformUtils initialization on task 0 only
  int ierr = 0;
  if ( ctxt_.onpe0() )
  {
    try
    {
      XMLPlatformUtils::Initialize();
    }

    catch (const XMLException& toCatch)
    {
      cout << "  <!-- Sample::readSample: Error during XML initialization :\n"
           << StrX(toCatch.getMessage()) << " -->" << endl;
      ierr = 1;
    }
    ctxt_.ibcast_send(1,1,&ierr,1);
  }
  else
  {
    ctxt_.ibcast_recv(1,1,&ierr,1,0,0);
  }
  //cout << ctxt_.mype() <<": SampleReader: ierr=" << ierr << endl;
  if ( ierr > 0 )
    throw SampleReaderException("error in XMLPlatformUtils::Initialize");
    
  // Determine if uri is a local file
  struct stat statbuf;
  bool read_from_file = !stat(uri.c_str(),&statbuf);
  
  // check for serial override
  read_from_file &= !serial;
  
  string xmlcontent;
  DoubleMatrix gfdata(ctxt_);
  if ( read_from_file )
  {
    const char* const filename = uri.c_str();
    if ( ctxt_.onpe0() )
      cout << " SampleReader: reading from file "
           << filename << " size: "
           << statbuf.st_size << endl;
         
    XMLGFPreprocessor xmlgfp;
    
    xmlgfp.process(filename,gfdata,xmlcontent);
    
    if ( ctxt_.onpe0() )
    {
      // cout << ctxt_.mype() << ": xmlcontent.size(): " << xmlcontent.size() 
      //      << endl;
      // cout << ctxt_.mype() << ": xmlcontent: " << xmlcontent << endl;
      memBufIS = new MemBufInputSource
        ( (const XMLByte*) &xmlcontent[0], xmlcontent.size(), "buf_id", false);
    }
  }
    
  bool read_wf = false;
  bool read_wfv = false;
  // initialize wavefunction_velocity in case it appears in the sample file
  
  if ( ctxt_.onpe0() )
  {
    parser = XMLReaderFactory::createXMLReader();
    if (valScheme == SAX2XMLReader::Val_Auto)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
        parser->setFeature(XMLUni::fgXercesDynamic, true);
    }

    if (valScheme == SAX2XMLReader::Val_Never)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, false);
    }

    if (valScheme == SAX2XMLReader::Val_Always)
    {
        parser->setFeature(XMLUni::fgSAX2CoreValidation, true);
        parser->setFeature(XMLUni::fgXercesDynamic, false);
    }

    parser->setFeature(XMLUni::fgSAX2CoreNameSpaces, doNamespaces);
    parser->setFeature(XMLUni::fgXercesSchema, doSchema);
    parser->setFeature(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);
    parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes, namespacePrefixes);

    int errorCount = 0;
    SampleHandler* s_handler = new SampleHandler(s,gfdata,wfvtmp);
    
    try
    {
      StructuredDocumentHandler handler(s_handler);
      parser->setContentHandler(&handler);
      parser->setErrorHandler(&handler);

#if 0
      //!!LocalFileInputSource localfile(XMLString::transcode(uri.c_str()));
      //!!parser->parse(localfile);
#endif

      cout << " Starting XML parsing" << endl;
      if ( read_from_file )
        parser->parse(*memBufIS);
      else
        parser->parse(uri.c_str());
      cout << " XML parsing done" << endl;
 
      errorCount = parser->getErrorCount();
    }

    catch (const XMLException& toCatch)
    {
        cout << "\nAn error occurred\n  Error: "
             << StrX(toCatch.getMessage())
             << "\n" << endl;
        //XMLPlatformUtils::Terminate();
        //delete parser;
        //throw;
    }

    catch (const SAXParseException& e)
    {
        cout << "\na SAXParseException occurred in file "
             << StrX(e.getSystemId())
             << ", line " << e.getLineNumber()
             << ", char " << e.getColumnNumber()
             << "\n  Message: " << StrX(e.getMessage()) << endl;
        //XMLPlatformUtils::Terminate();
        //delete parser;
        //throw;
    }
    
    read_wf = s_handler->read_wf;
    read_wfv = s_handler->read_wfv;
      
    delete s_handler;
    delete parser;
    XMLPlatformUtils::Terminate();
    
    // parsing of sample is complete, send end of sample tag to tasks > 0
    int tag = 0;
    ctxt_.ibcast_send(1,1,&tag,1);
    
  } // onpe0
  else
  {
    // tasks > 0
    // listen for Sample events
    // cout << ctxt_.mype() << ": listening..." << endl;
    
    bool done = false;
    int tag = -1;
    while ( !done )
    {
      ctxt_.ibcast_recv(1,1,&tag,1,0,0);
      // cout << ctxt_.mype() << ": received tag " << tag << endl;
      if ( tag == 0 )
        done = true;
      else if ( tag == 1 )
      {
        // Species
        string curr_name;
        ctxt_.string_bcast(curr_name,0);
        //cout << ctxt_.mype() << ": receiving species " << curr_name << endl;
        Species* sp = new Species(ctxt_,curr_name);
        SpeciesReader sp_reader(ctxt_);
        sp_reader.bcastSpecies(*sp);
        s.atoms.addSpecies(sp,curr_name);
      }        
      else if ( tag == 2 )
      {
        // Atom
        string curr_atom_name,curr_atom_species;
        D3vector curr_position, curr_velocity;
        ctxt_.string_bcast(curr_atom_name,0);
        ctxt_.string_bcast(curr_atom_species,0);
        double buf[3];
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        curr_position = D3vector(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        curr_velocity = D3vector(buf[0],buf[1],buf[2]);
        //cout << ctxt_.mype() << ": receiving atom " << curr_atom_name << endl;
        Atom* a = new Atom(curr_atom_name, curr_atom_species,
                           curr_position, curr_velocity);
        s.atoms.addAtom(a);
      }
      else if ( tag == 3 )
      {
        // Wavefunction
        read_wf = true;
        int nel,nspin,nempty;
        ctxt_.ibcast_recv(1,1,&nel,1,0,0);
        ctxt_.ibcast_recv(1,1,&nspin,1,0,0);
        ctxt_.ibcast_recv(1,1,&nempty,1,0,0);
        //cout << ctxt_.mype() << ": receiving wf nel=" << nel
        //     << " nspin=" << nspin << " nempty=" << nempty << endl;
        s.wf.set_nel(nel);
        s.wf.set_nspin(nspin);
        s.wf.set_nempty(nempty);
        
        // domain
        double buf[3];
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector a(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector b(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector c(buf[0],buf[1],buf[2]);
        UnitCell uc(a,b,c);
        
        // grid
        // receive only computed ecut
        double ecut;
        ctxt_.dbcast_recv(1,1,&ecut,1,0,0);
        
        s.wf.resize(uc,uc,ecut);
        
        //cout << ctxt_.mype() << ": wf resized, ecut=" << ecut << endl;
        
        //!! nkpoint fixed = 1
        const int nkpoint = 1;
        
        for ( int ispin = 0; ispin < nspin; ispin++ )
        {
          for ( int ikp = 0; ikp < nkpoint; ikp++ )
          {
            SlaterDet* sd = s.wf.sd(ispin,ikp);
            if ( sd != 0 )
            {
              // receive density_matrix
              vector<double> dmat_tmp(sd->nst());
              sd->context().dbcast_recv(sd->nst(),1,&dmat_tmp[0],1,0,0);
              sd->set_occ(dmat_tmp);
                
              if ( !read_from_file )
              {
                const Basis& basis = sd->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
                valarray<double> wftmpr(ft.np012loc());
                vector<complex<double> > wftmp(ft.np012loc());
                int size = -1;
                for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
                {
                  // read grid_function nloc fragment
                  //cout << sd->context();
                  //cout << sd->context().mype()
                  //     << ": wf receiving nloc=" << nloc << endl;
                  sd->context().irecv(1,1,&size,1,0,0);
                  //cout << sd->context().mype()
                  //     << ": received size=" << size << endl;
                  assert(size==ft.np012loc());
                  sd->context().drecv(size,1,&wftmpr[0],1,0,0);
                  //cout << sd->context().mype()
                  //     << ": grid_function nloc=" << nloc
                  //     << "received" << endl;
 
                  // copy to complex array
                  for ( int i = 0; i < size; i++ )
                  {
                    wftmp[i] = wftmpr[i];
                  }

                  ComplexMatrix& c = sd->c();
                  ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
                  //cout << sd->context().mype()
                  //     << ": grid_function read for nloc=" << nloc << endl;
                }
              }
            }
          }
        }
      }
      else if ( tag == 4 )
      {
        // cout << ctxt_.mype() << ": SampleReader received tag 4" << endl;
        read_wfv = true;
        // Wavefunction velocity
        int nel,nspin,nempty;
        ctxt_.ibcast_recv(1,1,&nel,1,0,0);
        ctxt_.ibcast_recv(1,1,&nspin,1,0,0);
        ctxt_.ibcast_recv(1,1,&nempty,1,0,0);
        //cout << ctxt_.mype() << ": receiving wf nel=" << nel
        //     << " nspin=" << nspin << " nempty=" << nempty << endl;
        wfvtmp.set_nel(nel);
        wfvtmp.set_nspin(nspin);
        wfvtmp.set_nempty(nempty);
        
        // domain
        double buf[3];
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector a(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector b(buf[0],buf[1],buf[2]);
        ctxt_.dbcast_recv(3,1,buf,1,0,0);
        D3vector c(buf[0],buf[1],buf[2]);
        UnitCell uc(a,b,c);
        
        // grid
        // receive only computed ecut
        double ecut;
        ctxt_.dbcast_recv(1,1,&ecut,1,0,0);
        
        wfvtmp.resize(uc,uc,ecut);
        
        //cout << ctxt_.mype() << ": wfv resized, ecut=" << ecut << endl;
        
        //!! nkpoint fixed = 1
        const int nkpoint = 1;
        
        for ( int ispin = 0; ispin < nspin; ispin++ )
        {
          for ( int ikp = 0; ikp < nkpoint; ikp++ )
          {
            SlaterDet* sd = wfvtmp.sd(ispin,ikp);
            if ( sd != 0 )
            {
              // receive density_matrix
              vector<double> dmat_tmp(sd->nst());
	      cout << " SampleReader: dmat.size()=" << dmat_tmp.size() << endl;
              sd->context().dbcast_recv(sd->nst(),1,&dmat_tmp[0],1,0,0);
              sd->set_occ(dmat_tmp);
                
              if ( !read_from_file )
              {
                const Basis& basis = sd->basis();
                FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
                valarray<double> wftmpr(ft.np012loc());
                vector<complex<double> > wftmp(ft.np012loc());
                int size = -1;
                for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
                {
                  // read grid_function nloc fragment
                  //cout << sd->context();
                  //cout << sd->context().mype()
                  //     << ": wfv receiving nloc=" << nloc << endl;
                  sd->context().irecv(1,1,&size,1,0,0);
                  //cout << sd->context().mype()
                  //     << ": received size=" << size << endl;
                  assert(size==ft.np012loc());
                  sd->context().drecv(size,1,&wftmpr[0],1,0,0);
                  //cout << sd->context().mype()
                  //     << ": grid_function nloc=" << nloc
                  //     << "received" << endl;
 
                  // copy to complex array
                  for ( int i = 0; i < size; i++ )
                  {
                    wftmp[i] = wftmpr[i];
                  }

                  ComplexMatrix& c = sd->c();
                  ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
                  //cout << sd->context().mype()
                  //     << ": grid_function read for nloc=" << nloc << endl;
                }
              }
            }
          }
        }
      }
    }
  } // if-else onpe0

  // This part is now executing on all tasks
  if ( read_from_file )
  {
    if ( read_wf )
    {
      // transfer data from the gfdata matrix to the SlaterDets
      //cout << ctxt_.mype() << ": mapping gfdata matrix..." 
      //     << endl;
      //cout << ctxt_.mype() << ": gfdata: (" << gfdata.m() << "x" << gfdata.n()
      //<< ") (" << gfdata.mb() << "x" << gfdata.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": gfdata.mloc()=" << gfdata.mloc()
      //<< " gfdata.nloc()=" << gfdata.nloc() << endl;
           
      //!! Next lines work for 1 kpoint, 1 spin only
      SlaterDet* sd = s.wf.sd(0,0);
      const Basis& basis = sd->basis();
      FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
      //cout << ctxt_.mype() << ": ft.np012loc()=" << ft.np012loc() << endl;
      //cout << ctxt_.mype() << ": ft.context(): " << ft.context();
      
      ComplexMatrix& c = sd->c();
      DoubleMatrix wftmpr(sd->context(),ft.np012(),sd->nst(),
        ft.np012loc(0),c.nb());
      //cout << ctxt_.mype() << ": wftmpr: (" << wftmpr.m() << "x" << wftmpr.n()
      //<< ") (" << wftmpr.mb() << "x" << wftmpr.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": wftmpr.mloc()=" << wftmpr.mloc()
      //<< " wftmpr.nloc()=" << wftmpr.nloc() << endl;

      vector<complex<double> > wftmp(ft.np012loc());      
      wftmpr.getsub(gfdata,ft.np012(),sd->nst(),0,0);
      
#if 0
      // Check orthogonality by computing overlap matrix
      DoubleMatrix smat(sd->context(),sd->nst(),sd->nst(),c.nb(),c.nb());
      smat.syrk('l','t',1.0,wftmpr,0.0);
      cout << smat;
#endif
      
      for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
      {
        // copy to complex array
        double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
        for ( int i = 0; i < ft.np012loc(); i++ )
          wftmp[i] = p[i];
        ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
      }
    } // if read_wf
    
    if ( read_wfv )
    {
      // transfer wfv data from the gfdata matrix to the SlaterDets
      //cout << ctxt_.mype() << ": mapping gfdata matrix..." 
      //     << endl;
      //cout << ctxt_.mype() << ": gfdata: (" << gfdata.m() << "x" << gfdata.n()
      //<< ") (" << gfdata.mb() << "x" << gfdata.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": gfdata.mloc()=" << gfdata.mloc()
      //<< " gfdata.nloc()=" << gfdata.nloc() << endl;
           
      //!! Next lines work for 1 kpoint, 1 spin only
      SlaterDet* sd = wfvtmp.sd(0,0);
      const Basis& basis = sd->basis();
      FourierTransform ft(basis,basis.np(0),basis.np(1),basis.np(2));
      //cout << ctxt_.mype() << ": ft.np012loc()=" << ft.np012loc() << endl;
      //cout << ctxt_.mype() << ": ft.context(): " << ft.context();
      
      ComplexMatrix& c = sd->c();
      DoubleMatrix wftmpr(sd->context(),ft.np012(),sd->nst(),
        ft.np012loc(0),c.nb());
      //cout << ctxt_.mype() << ": wftmpr: (" << wftmpr.m() << "x" << wftmpr.n()
      //<< ") (" << wftmpr.mb() << "x" << wftmpr.nb() << ") blocks" << endl;
      //cout << ctxt_.mype() << ": wftmpr.mloc()=" << wftmpr.mloc()
      //<< " wftmpr.nloc()=" << wftmpr.nloc() << endl;

      vector<complex<double> > wftmp(ft.np012loc());      
      wftmpr.getsub(gfdata,ft.np012(),sd->nst(),0,sd->nst());
      
      for ( int nloc = 0; nloc < sd->nstloc(); nloc++ )
      {
        // copy to complex array
        double* p = wftmpr.valptr(wftmpr.mloc()*nloc);
        for ( int i = 0; i < ft.np012loc(); i++ )
          wftmp[i] = p[i];
        ft.forward(&wftmp[0],c.valptr(c.mloc()*nloc));
      }
    }
  }
  
  // check if wavefunction_velocity element was read, if not, delete wfvtmp
  if ( s.wfv != 0 )
  {
    delete s.wfv;
    s.wfv = 0;
  }
  if ( read_wfv )
  {
    s.wfv = new Wavefunction(s.wf);
    *s.wfv = wfvtmp;
  }
}
