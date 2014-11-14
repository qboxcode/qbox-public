////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
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
// ElectricEnthalpy.C
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cmath>
#include <algorithm> // fill

#include "Timer.h"
#include "Context.h"
#include "Matrix.h"
#include "Sample.h"
#include "D3vector.h"
#include "ElectricEnthalpy.h"
#include "MLWFTransform.h"
#include "Wavefunction.h"
#include "D3vector.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "UnitCell.h"
#include "jade.h"
#include "blas.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::ElectricEnthalpy(Sample& s): s_(s), wf_(s.wf), dwf_(s.wf),
  sd_(*(s.wf.sd(0,0))), ctxt_(s.wf.sd(0,0)->context()),
  vctxt_(s.wf.sd(0,0)->basis().context())
{
  assert(wf_.nkp()==1);
  assert(wf_.nspin()==1);

  onpe0_ = ctxt_.onpe0();
  e_field_ = s.ctrl.e_field;
  compute_quadrupole_ = false;
  pol_type_ = s.ctrl.polarization_type;

  mlwft_ = new MLWFTransform(sd_);
  mlwft_->set_tol(1.e-10);

  vbasis_ = 0;
  smat_[0] = smat_[1] = smat_[2] = 0;
  rwf_[0] = rwf_[1] = rwf_[2] = 0;

  if ( (pol_type_ == "MLWF") || (pol_type_ == "MLWF_REF") )
  {
    // allocate real space wf arrays
    for ( int i = 0; i < 3; i++ )
      rwf_[i] = new Wavefunction(wf_);

    // Basis for real wavefunction
    vbasis_ = new Basis(vctxt_, D3vector(0,0,0));
    vbasis_->resize(wf_.cell(),wf_.refcell(),wf_.ecut()*4.0);
  }
  else if ( pol_type_ == "BERRY" )
  {
    // allocate complex Berry phase matrix
    int n = sd_.c().n();
    int nb = sd_.c().nb();
    for ( int i = 0; i < 3; i++ )
      smat_[i] = new ComplexMatrix(ctxt_,n,n,nb,nb);
  }
  else
  {
    cerr << "ElectricEnthalpy: invalid polarization type" << endl;
    ctxt_.abort(1);
  }

  niter_ = 1;

  if ( onpe0_ )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout.precision(8);
    cout << "<e_field> " << e_field_ << " </e_field>" << endl;
  }

  int nst = sd_.nst();
  mlwfc_.resize(nst);
  mlwfs_.resize(nst);
  correction_.resize(nst);
  correction_real_.resize(nst);
  quad_.resize(nst);
}

///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::~ElectricEnthalpy(void)
{
  delete mlwft_;
  delete vbasis_;
  for ( int i = 0; i < 3; i++ )
  {
    delete rwf_[i];
    delete smat_[i];
  }

  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<timing name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::update(void)
{
  // compute cos and sin matrices
  tmap["mlwf_update"].start();
  mlwft_->update();
  tmap["mlwf_update"].stop();
  vector<SlaterDet*> sdcos(3), sdsin(3);
  sdcos[0] = mlwft_->sdcosx();
  sdcos[1] = mlwft_->sdcosy();
  sdcos[2] = mlwft_->sdcosz();
  sdsin[0] = mlwft_->sdsinx();
  sdsin[1] = mlwft_->sdsiny();
  sdsin[2] = mlwft_->sdsinz();

  if ( pol_type_ == "MLWF" || pol_type_ == "MLWF_REF" )
  {
    tmap["mlwf_trans"].start();
    mlwft_->compute_transform();
    mlwft_->apply_transform(sd_);
    tmap["mlwf_trans"].stop();
    for ( int i = 0; i < sd_.nst(); i++ )
    {
      mlwfc_[i] = mlwft_->center(i);
      mlwfs_[i] = mlwft_->spread(i);
    }

    tmap["correction_real"].start();
    correction_real();
    tmap["correction_real"].stop();

    compute_polarization();

    // calculate gradient
    tmap["derivative_cor"].start();

    dwf_.clear();
    for ( int idir = 0; idir < 3; idir++ )
    {
      if ( e_field_[idir] != 0.0 )
      {
        // derivative of the electric enthalphy functional
        dwf_.clear();
        DoubleMatrix cc(rwf_[idir]->sd(0,0)->c());
        DoubleMatrix cp(dwf_.sd(0,0)->c());
        int size = cc.size();
        int ione = 1;
        double alpha = e_field_[idir];
        daxpy (&size, &alpha, cc.valptr(), &ione, cp.valptr(), &ione);
      }
    } //for idir
    tmap["derivative_cor"].stop();
  }
  else if ( pol_type_ == "BERRY" )
  {
    dwf_.clear();
    DoubleMatrix gradp(dwf_.sd(0,0)->c());

    for ( int idir = 0; idir < 3; idir++ )
    {
      if ( e_field_[idir] != 0.0 )
      {
        complex<double>* val = smat_[idir]->valptr();

        const double* re = mlwft_->a(idir*2)->cvalptr();
        const double* im = mlwft_->a(idir*2+1)->cvalptr();
        for ( int i = 0; i < smat_[idir]->size(); i++ )
          val[i] = complex<double>(re[i],im[i]);

  cout << "S[" << idir << "] = " << endl;
  cout << *smat_[idir] << endl;
        // invert S
        complex<double> z = smat_[idir]->inverse_det();
        gamma_[idir] = arg(z);

        int n = smat_[idir]->n();
        int nb = smat_[idir]->nb();

        // real and img part of S^{-1}
        DoubleMatrix s_real(ctxt_,n,n,nb,nb);
        DoubleMatrix s_img(ctxt_,n,n,nb,nb);
        DoubleMatrix sp(*smat_[idir]);

        int size = s_real.size();
        int ione = 1, itwo = 2;

        // copy real and imaginary parts of s to s_real and s_img
        dcopy(&size,sp.valptr(),&itwo,s_real.valptr(),&ione);
        dcopy(&size,sp.valptr()+1,&itwo,s_img.valptr(),&ione);

        // if ( ctxt_.onpe0() )cout << " S:\n";
        // cout << *smat_[idir];
        // if ( ctxt_.onpe0() )cout << " S_real:\n";
        // cout << s_real;
        // if ( ctxt_.onpe0() )cout << " S_img:\n";
        // cout << s_img;

        // proxy Matrix for cosx|psi> and sinx|psi>
        DoubleMatrix cosp(sdcos[idir]->c());
        DoubleMatrix sinp(sdsin[idir]->c());

        // alpha = E_i * L_i / 2pi
        //!! replace with reciprocal lattice vector
        const double fac = sd_.basis().cell().amat(idir*3+idir)/( 2.0*M_PI );
        double alpha = fac * e_field_[idir];

        gradp.gemm('n','n',alpha,cosp,s_img,1.0);
        gradp.gemm('n','n',alpha,sinp,s_real,1.0);
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_polarization_elec(void)
{
  polarization_elec_ = D3vector(0,0,0);
  polarization_elec_correction_ = D3vector(0,0,0);
  polarization_elec_correction_real_ = D3vector(0,0,0);
  for ( int i = 0; i < sd_.nst(); i++ )
  {
    polarization_elec_ -= 2.0 * mlwfc_[i];
    polarization_elec_correction_ -= 2.0 * correction_[i];
    polarization_elec_correction_real_ -= 2.0 * correction_real_[i];
  }
}
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_polarization_ion(void)
{
  polarization_ion_ = s_.atoms.dipole();
}
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_polarization(void)
{
  compute_polarization_ion();
  compute_polarization_elec();
  // total polarization
  polarization_ = polarization_ion_ + polarization_elec_;
}
///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::energy(Wavefunction& dwf, bool compute_hpsi)
{
  energy_ = - e_field_ * ( polarization_ + polarization_elec_correction_real_);
  if ( compute_hpsi )
  {
    // assert gamma-only and no spin
    assert(dwf.nkp()==1 && dwf.nspin()==1);
    dwf.sd(0,0)->c() += dwf_.sd(0,0)->c();
  }

  if ( onpe0_ )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout.precision(10);

    cout << " <polarization>\n";
    cout << "  <P_ion>    " << setw(18)
         << polarization_ion_  << " </P_ion>\n";
    cout << "  <P_elec>   " << setw(18)
         << polarization_elec_ << " </P_elec>\n";
    cout << "  <P_correct>" << setw(18)
         << polarization_elec_correction_ << " </P_correct>\n";
    cout << "  <P_cor_re> " << setw(18)
         << polarization_elec_correction_real_ << " </P_cor_re>\n";
    cout << "  <P_tot>    " << setw(18)
         << polarization_      << " </P_tot>\n";
    cout << "  <P_tot_cor>" << setw(18) << polarization_
            + polarization_elec_correction_ << " </P_tot_cor>\n";
    cout << "  <P_tot_re> " << setw(18) << polarization_
            + polarization_elec_correction_real_ << " </P_tot_re>\n";
    cout << " </polarization>\n";

    if ( compute_quadrupole_ )
    {
      //compute quadrupole associated with the mlwf center
      D3tensor q_mlwfc;
      D3tensor q_mlwfs;
      for ( int ist = 0; ist < sd_.nst(); ist++ )
      {
        D3vector& ctr = mlwfc_[ist];
        for (int j=0; j<3; j++)
        {
          for (int k = 0; k < 3; k++)
            q_mlwfc[j*3+k] -= 2.0 * ctr[j] * ctr[k];
        }

        q_mlwfs -= quad_[ist] * 2.0;
      }// for ist

      D3tensor q_ion = s_.atoms.quadrupole();
      D3tensor q_mlwf = q_mlwfc + q_mlwfs;
      //total primitive quadrupoles
      D3tensor q_tot = q_ion + q_mlwf;
      //traceless quadrupole
      D3tensor q_traceless = q_tot;
      q_traceless.traceless();

      cout << " <ionic_quadrupole> " << endl
           << q_ion
           << " </ionic_quadrupole>" << endl;
      cout << " <mlwfc_quadrupole> " << endl
           << q_mlwfc
           << " </mlwfc_quadrupole>" << endl;
      cout << " <mlwfs_quadrupole> " << endl
           << q_mlwfs
           << " </mlwfs_quadrupole>" << endl;
      cout << " <electronic_quadrupole> " << endl
           << q_mlwf
           << " </electronic_quadrupole>" << endl;
      cout << " <total_quadrupole> " << endl
           << q_tot
           << " </total_quadrupole>" << endl;
      cout << " <traceless_quadrupole> " << endl
           << q_traceless
           << " </traceless_quadrupole>" << endl;
      char uplo = 'u';
      D3vector eigval;
      D3tensor eigvec;
      q_traceless.syev(uplo, eigval, eigvec);
      cout << " <traceless_quadrupole_eigval> " << endl
           << eigval << endl
           << " </traceless_quadrupole_eigval>" << endl;
      cout << " <traceless_quadrupole_eigvec> " << endl
           << eigvec
           << " </traceless_quadrupole_eigvec>" << endl;

    } // if compute_quadrupole
  } // if onpe0
  return energy_;
}

///////////////////////////////////////////////////////////////////////////////
// Correction scheme by M. Stengel and N. Spaldin,
// Phys. Rev. B 73, 075121 (2006)
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::correction(void)
{
  const Basis& basis = sd_.basis();
  const Basis& vbasis = *vbasis_;
  int np0v = vbasis.np(0);
  int np1v = vbasis.np(1);
  int np2v = vbasis.np(2);
  const ComplexMatrix& c = sd_.c();
  DoubleMatrix cp(c);

  FourierTransform ft(basis,np0v,np1v,np2v);
  FourierTransform vft(vbasis,np0v,np1v,np2v);

  int np012loc = ft.np012loc();       //local rho grid size
  int nst = sd_.nst();                //total states
  int nloc = c.nloc();                //local states
  int mloc = c.mloc();                //wf(G) local size;
  int localsize= vbasis.localsize();  //rho(G) local size;


  correction_.resize(nst,D3vector(0,0,0));
  vector<double> ref(nst*3,0);

  for ( int in = 0; in < nloc; in++ )
  {
    //global state index
    int n = c.jglobal(in);

    vector<complex<double> > rhotmp(np012loc);
    vector<complex<double> > rhog(localsize);

    // wf in real space
    ft.backward(c.cvalptr(mloc*in),&rhotmp[0]);

    // rho in real space
    //#pragma omp parallel for
    for ( int i = 0; i < np012loc; i++ )
      rhotmp[i] = rhotmp[i] * rhotmp[i];

    // rho in Fourier space
    vft.forward(&rhotmp[0],&rhog[0]);

    double wb[3];
    for ( int idir = 0; idir < 3; idir++ )
    {
      wb[idir] = sd_.basis().cell().amat(idir*3+idir) / ( 2.0 * M_PI );
    }

    if (onpe0_) cout << "L"<< in << "  " << mlwfc_[n] << "   "
                     << wb[0] << "  " << wb[1] << "   "<< wb[2]<< endl;

    for ( int niter = 0; niter < niter_; niter++ )
    {
      // corrrection
      vector<double> dr(3,0);

      // correction for three Cartesian directions
      for ( int i = 0; i < localsize; i++ )
      {
        double gx[3];

        // determine if the G vector is on one and only one of the three axis;
        int naxis = 0, direc;
        for ( int idir = 0; idir < 3; idir++ )
        {
          gx[idir] = vbasis.gx( idir * localsize + i );

          if ( gx[idir] == 0.0 )
            naxis++;
          else
            direc = idir;
        }
        if ( naxis != 2 ) continue;

        int k = round ( gx[direc] * wb[direc] );
        complex<double> imag_one(0,1);

        //dr[direc] += -2.0 * volume * wb[direc] * real ( imag_one
        dr[direc] += -2.0 * wb[direc] * real ( imag_one * pow(-1.0,k)
                         / (double)k * rhog[i] * exp ( imag_one *
                          ( mlwfc_[n][direc] + ref[n*3+direc] )
                         * gx[direc] ) );
        /*
        if ( onpe0_ && i < 100 )
          cout << "N" << direc << "    " << k << "     "
               << gx[0] << "    " << gx[1] << "     " << gx[2] << "     "
               << dr[direc] << endl;
        */
      }// for i

      vctxt_.dsum(3,1,&dr[0],3);

      for ( int idir = 0; idir < 3; idir++ )
        ref[idir+3*n] += dr[idir];

      if ( onpe0_ ) cout << niter << ": " << ref[3*n] << "   "
                         << ref[3*n+1] << "   " << ref[3*n+2] << endl;

    }//niter

  }// for in


  ctxt_.dsum('r',nst*3,1,&ref[0],nst*3);

  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < nst; j++ )
      correction_[j][i] = ref[3*j+i];

}

////////////////////////////////////////////////////////////////////////////////
// Calculate correction in real space
// and derivatives of the correction
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::correction_real(void)
{
  const Basis& basis = sd_.basis();
  const Basis& vbasis = *vbasis_;
  int np0v = vbasis.np(0);
  int np1v = vbasis.np(1);
  int np2v = vbasis.np(2);
  const ComplexMatrix& c = sd_.c();
  DoubleMatrix cp(c);

  FourierTransform ft(basis,np0v,np1v,np2v);

  int np012v = ft.np012();
  int np012loc = ft.np012loc();
  int nst = sd_.nst();
  int nloc = c.nloc();
  int mloc = c.mloc();
  int ione = 1;
  // const double volume = vbasis.cell().volume();

  // store (x-x0)|psi> in rwfs
  rwf_[0]->clear();
  rwf_[1]->clear();
  rwf_[2]->clear();

  ComplexMatrix& cx = rwf_[0]->sd(0,0)->c();
  ComplexMatrix& cy = rwf_[1]->sd(0,0)->c();
  ComplexMatrix& cz = rwf_[2]->sd(0,0)->c();

  DoubleMatrix cpx(cx);
  DoubleMatrix cpy(cy);
  DoubleMatrix cpz(cz);

  // calculate refinements
  // ref is scaled by np012v
  vector<double> ref(nst*3);
  if ( compute_quadrupole_ ) ref.resize(nst*9);

  // cell size;
  const double ax = sd_.basis().cell().amat(0);
  const double ay = sd_.basis().cell().amat(4);
  const double az = sd_.basis().cell().amat(8);

  // half cell size;
  const double ax2 = ax / 2.0;
  const double ay2 = ay / 2.0;
  const double az2 = az / 2.0;

  // grid size;
  const double dx = ax / np0v;
  const double dy = ay / np1v;
  const double dz = az / np2v;

  /*
  for ( int i = 0; i < ctxt_.size(); i++ )
  {
    ctxt_.barrier();
    if ( i == ctxt_.myproc() )
      cout << i << " " << ctxt_.myrow() << "  " << ctxt_.mycol() << "  "
           << np012loc << endl;
  }
  */

  for ( int i = 0; i < nst; i++ )
    correction_real_[i] = D3vector(0,0,0);

  for ( int iter = 0; iter < niter_; iter++ )
  {
    fill(ref.begin(),ref.end(),0.0);

    for ( int in = 0; in < nloc; in++ )
    {
      int n = c.jglobal(in);
      //wavefunction in the real space
      vector<complex<double> > wftmp (np012loc,0);
      vector<complex<double> > xwftmp(np012loc,0);
      vector<complex<double> > ywftmp(np012loc,0);
      vector<complex<double> > zwftmp(np012loc,0);

      double* pref;
      if ( compute_quadrupole_ )
        pref = &ref[9*n];
      else
        pref = &ref[3*n];

      //real space wavefunction in wftmp
      tmap["ft"].start();
      ft.backward(c.cvalptr(mloc*in),&wftmp[0]);
      tmap["ft"].stop();

      tmap["real"].start();

      double x0 = mlwfc_[n][0] + correction_real_[n][0];
      double y0 = mlwfc_[n][1] + correction_real_[n][1];
      double z0 = mlwfc_[n][2] + correction_real_[n][2];
      //#pragma omp parallel for
      for ( int i = 0; i < np012loc; i++ )
      {
        double& pwf = *((double*)&wftmp[i]);
        double& pxwf = *((double*)&xwftmp[i]);
        double& pywf = *((double*)&ywftmp[i]);
        double& pzwf = *((double*)&zwftmp[i]);

        // do x;
        double x = dx * ( ft.i(i) ) - x0;
        if ( x > ax2 ) x -= ax;
        if ( x < -ax2 ) x += ax;

        pxwf = pwf * x;
        pref[0] += pwf * pxwf;

        // do y;
        double y = dy * ( ft.j(i) ) - y0;
        if ( y > ay2 ) y -= ay;
        if ( y < -ay2 ) y += ay;

        pywf = pwf * y;
        pref[1] += pwf * pywf;

        // do z;
        double z = dz * ( ft.k(i) ) - z0;
        if ( z > az2 ) z -= az;
        if ( z < -az2 ) z += az;

        pzwf = pwf * z;
        pref[2] += pwf * pzwf;

        if ( compute_quadrupole_ )
        {
          pref[3] += pxwf * pxwf;
          pref[4] += pywf * pywf;
          pref[5] += pzwf * pzwf;
          pref[6] += pxwf * pywf;
          pref[7] += pywf * pzwf;
          pref[8] += pzwf * pxwf;
        }
        /*
        if (onpe0_&& i< 100 )
          cout << i << "  " << pwf << "   " << pxwf << "  "
               << pref[0] << "  " << pref[1] << "  " << pref[2] << endl;
        */
      } //for i

      /*
      cout.precision(10);
      if ( onpe0_ ) cout << "JJ" << n << "  " << iter << "  "
                         << pref[0] << "  " << pref[1] << "  " << pref[2] << "  "
                         << pref[3] << "  " << pref[4] << "  " << pref[5] << "  "
                         << endl;

      */
      tmap["real"].stop();

      // ft to get xwf in the reciprocal space at the last iteration
      if ( iter == niter_ - 1 )
      {
        tmap["ft"].start();
        if ( e_field_[0] != 0.0 ) ft.forward(&xwftmp[0],cx.valptr(mloc*in));
        if ( e_field_[1] != 0.0 ) ft.forward(&ywftmp[0],cy.valptr(mloc*in));
        if ( e_field_[2] != 0.0 ) ft.forward(&zwftmp[0],cz.valptr(mloc*in));
        tmap["ft"].stop();
      } // if

    } //for in

    ctxt_.barrier();
    tmap["dsum"].start();
    if ( compute_quadrupole_ )
      ctxt_.dsum(9*nst,1,&ref[0],9*nst);
    else
      ctxt_.dsum(3*nst,1,&ref[0],3*nst);
    tmap["dsum"].stop();


    tmap["real"].start();

    if ( compute_quadrupole_ )
    {
      #pragma omp parallel for
      for ( int ist = 0; ist < nst; ist++ )
      {
        D3vector& pcor = correction_real_[ist];
        D3tensor& pquad = quad_[ist];

        pcor[0] += ref[ist*9]/np012v;
        pcor[1] += ref[ist*9+1]/np012v;
        pcor[2] += ref[ist*9+2]/np012v;
        pquad.setdiag ( 0, ref[ist*9+3]/np012v - pcor[0] * pcor[0] );
        pquad.setdiag ( 1, ref[ist*9+4]/np012v - pcor[1] * pcor[1] );
        pquad.setdiag ( 2, ref[ist*9+5]/np012v - pcor[2] * pcor[2] );
        pquad.setoffdiag ( 0, ref[ist*9+6]/np012v - pcor[0] * pcor[1] );
        pquad.setoffdiag ( 1, ref[ist*9+7]/np012v - pcor[1] * pcor[2] );
        pquad.setoffdiag ( 2, ref[ist*9+8]/np012v - pcor[2] * pcor[0] );

      }
    }
    else
    {
      #pragma omp parallel for
      for ( int ist = 0; ist < nst; ist++ )
      {
        D3vector& pcor = correction_real_[ist];

        pcor[0] += ref[ist*3]/np012v;
        pcor[1] += ref[ist*3+1]/np012v;
        pcor[2] += ref[ist*3+2]/np012v;
      }
    }
    tmap["real"].stop();

  } //for iter
}
