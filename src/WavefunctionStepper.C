////////////////////////////////////////////////////////////////////////////////
//
// WavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: WavefunctionStepper.C,v 1.1 2003-11-21 20:01:06 fgygi Exp $

#include "WavefunctionStepper.h"
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
using namespace std;

#if 0
////////////////////////////////////////////////////////////////////////////////
void WavefunctionStepper::diag(Wavefunction& dwf)
{
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      if ( wf.sd(ispin,ikp) != 0 )
      {
        if ( wf.sdcontext(ispin,ikp)->active() )
        {
          if ( s_.ctrl.wf_diag == "T" )
          {
          tmap["eigval"].start();
          // compute eigenvalues
          if ( wf.sd(ispin,ikp)->basis().real() )
          {
            // proxy real matrices c, cp
            DoubleMatrix c(wf.sd(ispin,ikp)->c());
            DoubleMatrix cp(dwf.sd(ispin,ikp)->c());
 
            DoubleMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
            // factor 2.0 in next line: G and -G
            h.gemm('t','n',2.0,c,cp,0.0);
            // rank-1 update correction
            h.ger(-1.0,c,0,cp,0);
 
            //cout << " Hamiltonian at k = " << wf.sd(ispin,ikp)->kpoint()
            //     << endl;
            //cout << h;
 
            valarray<double> w(h.m());
            h.syev('l',w);
            if ( s_.wf.context().onpe0() )
            {
              const double eVolt = 2.0 * 13.6058;
              cout <<    "  <eigenvalues spin=\"" << ispin
                   << "\" kpoint=\"" << wf.sd(ispin,ikp)->kpoint()
                   << "\" n=\"" << h.m() << "\">" << endl;
              for ( int i = 0; i < h.m(); i++ )
              {
                cout << setw(10) << setprecision(5) << w[i]*eVolt;
                if ( i%5 == 4 ) cout << endl;
              }
              if ( h.m()%5 != 0 ) cout << endl;
              cout << "  </eigenvalues>" << endl;
            }
          }
          else
          {
            // complex case not implemented
            assert(false);
            #if 0
            ComplexMatrix& c(wf.sd[ikp]->c());
            ComplexMatrix& cp(dwf.sd[ikp]->c());
 
            ComplexMatrix h(c.context(),c.n(),c.n(),c.nb(),c.nb());
 
            h.gemm('c','n',1.0,c,cp,0.0);
 
            //cout << " Hamiltonian at k = " << wf.sd[ikp]->kpoint() << endl;
            //cout << h;
 
            valarray<double> w(h.m());

            h.heev('l',w);
            cout << " Eigenvalues at k = " << wf.sd[ikp]->kpoint() << endl;
            const double eVolt = 2.0 * 13.6058;
            for ( int i = 0; i < h.m(); i++ )
            {
              cout << "%" << setw(3) << ikp
                   << setw(10) << setprecision(5) << w[i]*eVolt << endl;;
            }
            #endif
          }
          tmap["eigval"].stop();
 
          } // wfdiag T
        }
      }
    }
  }
}

#endif
