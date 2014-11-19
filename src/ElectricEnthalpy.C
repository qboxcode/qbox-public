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

  if ( s.ctrl.polarization_type == "BERRY" )
    pol_type_ = berry;
  else if ( s.ctrl.polarization_type == "MLWF" )
    pol_type_ = mlwf;
  else if ( s.ctrl.polarization_type == "MLWF_REF" )
    pol_type_ = mlwf_ref;
  else
  {
    cerr << "ElectricEnthalpy: invalid polarization type" << endl;
    ctxt_.abort(1);
  }

  mlwft_ = new MLWFTransform(sd_);
  mlwft_->set_tol(1.e-10);

  vbasis_ = 0;
  smat_[0] = smat_[1] = smat_[2] = 0;
  rwf_[0] = rwf_[1] = rwf_[2] = 0;
  int nst = sd_.nst();

  if ( pol_type_ == mlwf_ref )
  {
    // allocate real space wf arrays for MLWF refinement
    for ( int i = 0; i < 3; i++ )
      rwf_[i] = new Wavefunction(wf_);

    // Basis for real wavefunction
    vbasis_ = new Basis(vctxt_, D3vector(0,0,0));
    vbasis_->resize(wf_.cell(),wf_.refcell(),wf_.ecut()*4.0);
    correction_.resize(nst);
  }
  else if ( pol_type_ == berry )
  {
    // allocate complex Berry phase matrix
    int n = sd_.c().n();
    int nb = sd_.c().nb();
    for ( int i = 0; i < 3; i++ )
      smat_[i] = new ComplexMatrix(ctxt_,n,n,nb,nb);
  }

  if ( onpe0_ )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout.precision(8);
    cout << "<e_field> " << e_field_ << " </e_field>" << endl;
  }

  mlwfc_.resize(nst);
  mlwfs_.resize(nst);
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

  polarization_ion_ = s_.atoms.dipole();
  polarization_elec_ = D3vector(0,0,0);

  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref )
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

    if ( pol_type_ == mlwf_ref )
    {
      tmap["correction"].start();
      compute_correction();
      tmap["correction"].stop();
    }

    // calculate gradient
    dwf_.clear();
    for ( int idir = 0; idir < 3; idir++ )
    {
      if ( e_field_[idir] != 0.0 )
      {
        // MLWF part
        if ( pol_type_ == mlwf )
        {
          const double nst = sd_.nst();
          std::vector<double> adiag_inv_real(nst,0),adiag_inv_imag(nst,0);
          for ( int ist = 0; ist < nst; ist ++ )
          {
            const std::complex<double>
              adiag( mlwft_->adiag(idir*2,ist),mlwft_->adiag(idir*2+1,ist) );

            adiag_inv_real[ist] = real( std::complex<double>(1,0) / adiag );
            adiag_inv_imag[ist] = imag( std::complex<double>(1,0) / adiag );

          }

          DoubleMatrix ccos(sdcos[idir]->c());
          DoubleMatrix csin(sdsin[idir]->c());
          DoubleMatrix cp(dwf_.sd(0,0)->c());

          int nloc = cp.nloc();
          int mloc = cp.mloc();
          int ione = 1;

          const double fac = sd_.basis().cell().amat(idir*3+idir)
                             * e_field_[idir] / ( 2.0 * M_PI );

          for (int in = 0; in < nloc; in++)
          {
            int ist = cp.jglobal(in);
            double fac1 = adiag_inv_real[ist] * fac;
            double fac2 = adiag_inv_imag[ist] * fac;

            double *ptr1 = &cp[in*mloc],
                   *ptrcos = &ccos[in*mloc],
                   *ptrsin = &csin[in*mloc];

            daxpy(&mloc, &fac2, ptrcos, &ione, ptr1, &ione);
            daxpy(&mloc, &fac1, ptrsin, &ione, ptr1, &ione);

          }
        }
        else if ( pol_type_ == mlwf_ref )
        {
          // MLWF_REF part: real-space correction
          DoubleMatrix cc(rwf_[idir]->sd(0,0)->c());
          DoubleMatrix cp(dwf_.sd(0,0)->c());

          int size = cc.size();
          double alpha = e_field_[idir];
          int ione = 1;
          daxpy (&size, &alpha, cc.valptr(), &ione, cp.valptr(), &ione);
        } // if pol_type_
      } // if e_field_[idir]
    } // for idir

    for ( int i = 0; i < sd_.nst(); i++ )
    {
      polarization_elec_ -= 2.0 * mlwfc_[i];
      if ( pol_type_ == mlwf_ref )
        polarization_elec_ -= 2.0 * correction_[i];
    }
  }
  else if ( pol_type_ == berry )
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

        // invert S and compute determinant
        complex<double> z = smat_[idir]->inverse_det();
        double gamma = arg(z);

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

        // assume occupation number of 2.0
        polarization_elec_[idir] = - 2.0 * fac * gamma;

      } // if e_field_[idir]
    }
  }

  polarization_ = polarization_ion_ + polarization_elec_;
}

///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::energy(Wavefunction& dwf, bool compute_hpsi)
{
  energy_ = - e_field_ * polarization_;
  if ( compute_hpsi )
  {
    // assert gamma-only and no spin
    assert(dwf.nkp()==1 && dwf.nspin()==1);
    dwf.sd(0,0)->c() += dwf_.sd(0,0)->c();
  }
  return energy_;
}

///////////////////////////////////////////////////////////////////////////////
// Correction scheme by M. Stengel and N. Spaldin,
// Phys. Rev. B 73, 075121 (2006)
// Calculate correction in real space and derivatives of the correction
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_correction(void)
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

  for ( int i = 0; i < nst; i++ )
    correction_[i] = D3vector(0,0,0);

  const int niter_ = 1;
  for ( int iter = 0; iter < niter_; iter++ )
  {
    fill(ref.begin(),ref.end(),0.0);

    for ( int in = 0; in < nloc; in++ )
    {
      int n = c.jglobal(in);
      // wavefunction in real space
      vector<complex<double> > wftmp (np012loc,0);
      vector<complex<double> > xwftmp(np012loc,0);
      vector<complex<double> > ywftmp(np012loc,0);
      vector<complex<double> > zwftmp(np012loc,0);

      double* pref;
      if ( compute_quadrupole_ )
        pref = &ref[9*n];
      else
        pref = &ref[3*n];

      // real space wavefunction in wftmp
      tmap["ft"].start();
      ft.backward(c.cvalptr(mloc*in),&wftmp[0]);
      tmap["ft"].stop();

      tmap["real"].start();
      double x0 = mlwfc_[n][0] + correction_[n][0];
      double y0 = mlwfc_[n][1] + correction_[n][1];
      double z0 = mlwfc_[n][2] + correction_[n][2];
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
      } // for i
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
        D3vector& pcor = correction_[ist];
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
      for ( int ist = 0; ist < nst; ist++ )
      {
        D3vector& pcor = correction_[ist];

        pcor[0] += ref[ist*3]/np012v;
        pcor[1] += ref[ist*3+1]/np012v;
        pcor[2] += ref[ist*3+2]/np012v;
      }
    }
    tmap["real"].stop();
  } // for iter
}

////////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::print(ostream& os) const
{
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  // print MLWF centers if pol_type_ == MLWF or MLWF_REF
  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref )
  {
    int nst = sd_.nst();
    os << " <mlwf_set size=\"" << nst << "\">" << endl;
    os << setprecision(8);
    for ( int i = 0; i < nst; i++ )
    {
      os << " <mlwf center=\"" << setprecision(8)
         << setw(12) << mlwfc_[i].x
         << setw(12) << mlwfc_[i].y
         << setw(12) << mlwfc_[i].z
         << " \" spread=\" " << mlwfs_[i] << " \"/>" << endl;
      if ( pol_type_ == mlwf_ref )
      {
        os << " <mlwf_ref center=\"" << setprecision(8)
           << setw(12) << mlwfc_[i].x + correction_[i].x
           << setw(12) << mlwfc_[i].y + correction_[i].y
           << setw(12) << mlwfc_[i].z + correction_[i].z;
        if ( compute_quadrupole_ )
        {
          // add spread attribute
          os << " \n\" spread=\" " << sqrt(quad_[i].trace()) << " \"";
        }
        os << "/>" << endl;

        if ( compute_quadrupole_ )
          os << "    <quad>"
             << setw(12) << quad_[i][0]
             << setw(12) << quad_[i][4]
             << setw(12) << quad_[i][8]
             << setw(12) << quad_[i][1]
             << setw(12) << quad_[i][2]
             << setw(12) << quad_[i][5]
             << " </quad>" << endl;
      }
    }
    os << " </mlwf_set>" << endl;

    //compute quadrupole associated with the mlwf center
    if ( compute_quadrupole_ )
    {
      D3tensor q_mlwfc;
      D3tensor q_mlwfs;
      for ( int ist = 0; ist < sd_.nst(); ist++ )
      {
        D3vector ctr = mlwfc_[ist];
        for (int j=0; j<3; j++)
        {
          for (int k = 0; k < 3; k++)
            q_mlwfc[j*3+k] -= 2.0 * ctr[j] * ctr[k];
        }
        q_mlwfs -= quad_[ist] * 2.0;
      } // for ist

      D3tensor q_ion = s_.atoms.quadrupole();
      D3tensor q_mlwf = q_mlwfc + q_mlwfs;
      //total primitive quadrupoles
      D3tensor q_tot = q_ion + q_mlwf;
      //traceless quadrupole
      D3tensor q_traceless = q_tot;
      q_traceless.traceless();

      os << " <ionic_quadrupole> " << endl
           << q_ion
           << " </ionic_quadrupole>" << endl;
      os << " <mlwfc_quadrupole> " << endl
           << q_mlwfc
           << " </mlwfc_quadrupole>" << endl;
      os << " <mlwfs_quadrupole> " << endl
           << q_mlwfs
           << " </mlwfs_quadrupole>" << endl;
      os << " <electronic_quadrupole> " << endl
           << q_mlwf
           << " </electronic_quadrupole>" << endl;
      os << " <total_quadrupole> " << endl
           << q_tot
           << " </total_quadrupole>" << endl;
      os << " <traceless_quadrupole> " << endl
           << q_traceless
           << " </traceless_quadrupole>" << endl;
      char uplo = 'u';
      D3vector eigval;
      D3tensor eigvec;
      q_traceless.syev(uplo, eigval, eigvec);
      os << " <traceless_quadrupole_eigval> " << endl
           << eigval << endl
           << " </traceless_quadrupole_eigval>" << endl;
      os << " <traceless_quadrupole_eigvec> " << endl
           << eigvec
           << " </traceless_quadrupole_eigvec>" << endl;
    }
  }

  // print polarization
  os.precision(10);
  os << "  <polarization>\n";
  os << "    <P_ion>  " << setw(16) << polarization_ion_.x
                        << setw(16) << polarization_ion_.y
                        << setw(16) << polarization_ion_.z << " </P_ion>\n";
  os << "    <P_elec> " << setw(16) << polarization_elec_.x
                        << setw(16) << polarization_elec_.y
                        << setw(16) << polarization_elec_.z << " </P_elec>\n";
  os << "    <P_tot>  " << setw(16) << polarization_.x
                        << setw(16) << polarization_.y
                        << setw(16) << polarization_.z << " </P_tot>\n";
  os << "  </polarization>\n";
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const ElectricEnthalpy& e )
{
  e.print(os);
  return os;
}
