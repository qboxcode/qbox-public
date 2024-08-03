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
// ElectricEnthalpy.cpp
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
using namespace std;

///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::vsst(double x) const
{
  // smooth sawtooth periodic potential function
  // x in [-1/2, 1/2]
  // The slope of vsst is 1 at x=0
  //
  // The function vsst approximates the identity function x->x
  // in the interval [-1/2+xcut, 1/2-xcut]
  const double xcut = 0.05;
  const double xcut2 = xcut*xcut;
  // The function vsst is well represented by a
  // discrete Fourier transform of length np = 2*ng
  // Some aliasing error will arise if np < 2*ng
  // The error is determined by the product xcut*ng
  const int ng = 32;
  double v = 0.0;
  for ( int ig = 1; ig < ng; ig++ )
  {
    const double g = 2 * M_PI * ig;
    const double arg = -0.25 * g * g * xcut2;
    // next line: factor sgn to shift origin by 0.5
    const int sgn = 1 - 2*(ig%2);
    const double c = -2.0 * sgn * exp ( arg ) / g;
    v += c * sin(x*g);
  }
  return v;
}

///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::ElectricEnthalpy(Sample& s): s_(s), wf_(s.wf),
  sd_(*(s.wf.sd(0,0))), ctxt_(s.wf.sd(0,0)->context()),
  basis_(s.wf.sd(0,0)->basis())
{
  onpe0_ = MPIdata::onpe0();
  e_field_ = s.ctrl.e_field;
  finite_field_ = norm2(e_field_) != 0.0;
  compute_quadrupole_ = false;

  if ( s.ctrl.polarization == "OFF" )
    pol_type_ = off;
  else if ( s.ctrl.polarization == "BERRY" )
    pol_type_ = berry;
  else if ( s.ctrl.polarization == "MLWF" )
    pol_type_ = mlwf;
  else if ( s.ctrl.polarization == "MLWF_REF" )
    pol_type_ = mlwf_ref;
  else if ( s.ctrl.polarization == "MLWF_REF_Q" )
  {
    pol_type_ = mlwf_ref_q;
    compute_quadrupole_ = true;
  }
  else
  {
    cerr << "ElectricEnthalpy: invalid polarization option" << endl;
    ctxt_.abort(1);
  }

  // do not allocate further objects if polarization is off
  if ( pol_type_ == off ) return;

  if ( wf_.nkp() != 1 )
    throw invalid_argument("ElectricEnthalpy: nkp != 1");
  if ( wf_.nspin() != 1 )
    throw invalid_argument("ElectricEnthalpy: nspin != 1");
  // check that there are no fractionally occupied states
  // next line: (3-nspin) = 2 if nspin==1 and 1 if nspin==2
  if ( wf_.nel() != ( 2 * wf_.nst() ) )
  {
    throw invalid_argument("ElectricEnthalpy: fractionally occupied"
      " or empty states");
  }

  dwf_ = new Wavefunction(s.wf);
  mlwft_ = new MLWFTransform(sd_);
  mlwft_->set_tol(1.e-10);

  smat_[0] = smat_[1] = smat_[2] = 0;
  rwf_[0] = rwf_[1] = rwf_[2] = 0;
  int nst = sd_.nst();

  if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
  {
    // allocate real space wf arrays for MLWF refinement
    for ( int idir = 0; idir < 3; idir++ )
      rwf_[idir] = new Wavefunction(wf_);
    correction_.resize(nst);
  }
  else if ( pol_type_ == berry )
  {
    // allocate complex Berry phase matrix
    int n = sd_.c().n();
    int nb = sd_.c().nb();
    for ( int idir = 0; idir < 3; idir++ )
      smat_[idir] = new ComplexMatrix(ctxt_,n,n,nb,nb);
  }

  if ( (pol_type_ != off) && onpe0_ )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout.precision(8);
    cout << "<e_field> " << e_field_ << " </e_field>" << endl;
  }

  mlwfc_.resize(nst);
  mlwfs_.resize(nst);
  quad_.resize(nst);

  tmap["mlwf_update"].reset();
  tmap["mlwf_trans"].reset();
  tmap["correction"].reset();
  tmap["ft"].reset();
  tmap["real"].reset();
  tmap["ft"].reset();
  tmap["dsum"].reset();
  tmap["real"].reset();
}

///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::~ElectricEnthalpy(void)
{
  if ( pol_type_ == off ) return;

  delete dwf_;
  delete mlwft_;
  for ( int idir = 0; idir < 3; idir++ )
  {
    delete rwf_[idir];
    delete smat_[idir];
  }

  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = i->second.real();
    double tmin, tmax;
    MPI_Reduce(&time,&tmin,1,MPI_DOUBLE,MPI_MIN,0,MPIdata::comm());
    MPI_Reduce(&time,&tmax,1,MPI_DOUBLE,MPI_MAX,0,MPIdata::comm());
    if ( MPIdata::onpe0() && (tmax > 0.0) )
    {
      string s = "name=\"" + (*i).first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::update(void)
{
  if ( pol_type_ == off ) return;

  const UnitCell& cell = sd_.basis().cell();
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

  dipole_ion_ = s_.atoms.dipole();
  dipole_el_ = D3vector(0,0,0);

  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
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

    if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
    {
      tmap["correction"].start();
      compute_correction();
      tmap["correction"].stop();
    }

    for ( int i = 0; i < sd_.nst(); i++ )
    {
      dipole_el_ -= 2.0 * mlwfc_[i];
      if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
        dipole_el_ -= 2.0 * correction_[i];
    }

    // compute gradient
    if ( finite_field_ )
    {
      dwf_->clear();
      for ( int idir = 0; idir < 3; idir++ )
      {
        // MLWF part
        if ( pol_type_ == mlwf )
        {
          const double fac = ( e_field_ * cell.a(idir) ) / ( 2.0 * M_PI );

          if ( fac != 0.0 )
          {
            const double nst = sd_.nst();
            std::vector<double> z_inv_real(nst,0),z_inv_imag(nst,0);
            for ( int ist = 0; ist < nst; ist ++ )
            {
              // c,s = B^T * adiag(idir,ist)
              const double itwopi = 1.0 / ( 2.0 * M_PI );
              const double *bmat = cell.bmat();
              double c,s;
              if ( idir == 0 )
              {
                c = itwopi * ( bmat[0] * mlwft_->adiag(0,ist) +
                               bmat[1] * mlwft_->adiag(2,ist) +
                               bmat[2] * mlwft_->adiag(4,ist) );
                s = itwopi * ( bmat[0] * mlwft_->adiag(1,ist) +
                               bmat[1] * mlwft_->adiag(3,ist) +
                               bmat[2] * mlwft_->adiag(5,ist) );
              }
              else if ( idir == 1 )
              {
                c = itwopi * ( bmat[3] * mlwft_->adiag(0,ist) +
                               bmat[4] * mlwft_->adiag(2,ist) +
                               bmat[5] * mlwft_->adiag(4,ist) );
                s = itwopi * ( bmat[3] * mlwft_->adiag(1,ist) +
                               bmat[4] * mlwft_->adiag(3,ist) +
                               bmat[5] * mlwft_->adiag(5,ist) );
              }
              else
              {
                c = itwopi * ( bmat[6] * mlwft_->adiag(0,ist) +
                               bmat[7] * mlwft_->adiag(2,ist) +
                               bmat[8] * mlwft_->adiag(4,ist) );
                s = itwopi * ( bmat[6] * mlwft_->adiag(1,ist) +
                               bmat[7] * mlwft_->adiag(3,ist) +
                               bmat[8] * mlwft_->adiag(5,ist) );
              }
              const std::complex<double> z(c,s);
              const std::complex<double> z_inv = std::complex<double>(1,0)/z;
              z_inv_real[ist] = real( z_inv );
              z_inv_imag[ist] = imag( z_inv );
            }

            DoubleMatrix ccos(sdcos[idir]->c());
            DoubleMatrix csin(sdsin[idir]->c());
            DoubleMatrix cp(dwf_->sd(0,0)->c());

            int nloc = cp.nloc();
            int mloc = cp.mloc();
            int ione = 1;


            for (int in = 0; in < nloc; in++)
            {
              int ist = cp.jglobal(in);
              double fac1 = z_inv_real[ist] * fac;
              double fac2 = z_inv_imag[ist] * fac;

              double *ptr1 = &cp[in*mloc],
                     *ptrcos = &ccos[in*mloc],
                     *ptrsin = &csin[in*mloc];

              daxpy(&mloc, &fac2, ptrcos, &ione, ptr1, &ione);
              daxpy(&mloc, &fac1, ptrsin, &ione, ptr1, &ione);
            }
          } // if fac
        }
        else if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
        {
          // MLWF_REF part: real-space correction
          if ( e_field_[idir] != 0.0 )
          {
            DoubleMatrix cc(rwf_[idir]->sd(0,0)->c());
            DoubleMatrix cp(dwf_->sd(0,0)->c());

            int size = cc.size();
            double alpha = e_field_[idir];
            int ione = 1;
            daxpy (&size, &alpha, cc.valptr(), &ione, cp.valptr(), &ione);
          }
        } // if pol_type_
      } // for idir
    } // if finite_field_
  }
  else if ( pol_type_ == berry )
  {
    // proxy matrix
    DoubleMatrix gradp(dwf_->sd(0,0)->c());
    if ( finite_field_ )
      dwf_->clear();

    for ( int idir = 0; idir < 3; idir++ )
    {
      complex<double>* val = smat_[idir]->valptr();

      const double* re = mlwft_->b(idir*2)->cvalptr();
      const double* im = mlwft_->b(idir*2+1)->cvalptr();
      for ( int i = 0; i < smat_[idir]->size(); i++ )
        val[i] = complex<double>(re[i],im[i]);

      // compute determinant of S
      ComplexMatrix& sm = *smat_[idir];
      const Context& ctxt = sm.context();

      // perform LU decomposition of S
      valarray<int> ipiv;
      sm.lu(ipiv);

      int n = sm.n();
      int nb = sm.nb();
      // note: m == n, mb == nb

      // compute determinant from diagonal values of U  (stored in diag of S)
      // get full diagonal of the matrix in array diag
      valarray<complex<double> > diag(n);
      for ( int ii = 0; ii < n; ii++ )
      {
        int iii = sm.l(ii) * nb + sm.x(ii);
        int jjj = sm.m(ii) * nb + sm.y(ii);
        if ( sm.pr(ii) == ctxt.myrow() &&
             sm.pc(ii) == ctxt.mycol() )
          diag[ii] = val[iii+sm.mloc()*jjj];
      }
      ctxt.dsum(2*n,1,(double*)&diag[0],2*n);

      // sum the argument of diagonal elements
      double sumarg = 0.0;
      for ( int ii = 0; ii < n; ii++ )
        sumarg += arg(diag[ii]);

      const int sign = sm.signature(ipiv);
      if ( sign == -1 )
        sumarg += M_PI;

      // assume occupation number of 2.0
      // sumarg contains Im ln S(idir)
      // contribution to dipole[i] = sum_idir -2.0 * a(i,idir) * S(idir)/(2*pi)

      dipole_el_ -= cell.a(idir) * sumarg / M_PI;

      if ( finite_field_ )
      {
        // alpha = ( e_field dot a_idir ) / ( 2 * pi )
        double alpha = e_field_ * cell.a(idir) / ( 2.0 * M_PI );

        if ( alpha != 0.0 )
        {

          // compute inverse of smat
          sm.inverse_from_lu(ipiv);

          // real and img part of S^{-1}
          DoubleMatrix s_real(ctxt_,n,n,nb,nb);
          DoubleMatrix s_img(ctxt_,n,n,nb,nb);
          DoubleMatrix sp(sm);

          int size = s_real.size();
          int ione = 1, itwo = 2;

          // copy real and imaginary parts of s to s_real and s_img
          dcopy(&size,sp.valptr(),&itwo,s_real.valptr(),&ione);
          dcopy(&size,sp.valptr()+1,&itwo,s_img.valptr(),&ione);

          // proxy Matrix for cosx|psi> and sinx|psi>
          DoubleMatrix cosp(sdcos[idir]->c());
          DoubleMatrix sinp(sdsin[idir]->c());

          gradp.gemm('n','n',alpha,cosp,s_img,1.0);
          gradp.gemm('n','n',alpha,sinp,s_real,1.0);
        }
      }
    } // for idir
  }

  dipole_total_ = dipole_ion_ + dipole_el_;
}

///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::enthalpy(Wavefunction& dwf, bool compute_hpsi)
{
  // return zero if polarization is off, or field is zero
  if ( pol_type_ == off || !finite_field_ )
    return 0.0;

  enthalpy_ = - e_field_ * dipole_total_;
  if ( compute_hpsi )
  {
    // assert gamma-only and no spin
    assert(dwf.nkp()==1 && dwf.nspin()==1);
    dwf.sd(0,0)->c() += dwf_->sd(0,0)->c();
  }
  return enthalpy_;
}

///////////////////////////////////////////////////////////////////////////////
// Correction scheme by M. Stengel and N. Spaldin,
// Phys. Rev. B 73, 075121 (2006)
// Calculate correction in real space and derivatives of the correction
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_correction(void)
{
  int np0v = basis_.np(0);
  int np1v = basis_.np(1);
  int np2v = basis_.np(2);
  const ComplexMatrix& c = sd_.c();
  DoubleMatrix cp(c);

  FourierTransform ft(basis_,np0v,np1v,np2v);

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

  // calculate refinements
  // ref is scaled by np012v
  vector<double> ref(nst*3);
  if ( compute_quadrupole_ ) ref.resize(nst*9);

  // cell size;
  const UnitCell& cell = sd_.basis().cell();
  const D3vector a0 = cell.a(0);
  const D3vector a1 = cell.a(1);
  const D3vector a2 = cell.a(2);

  // cell vector lengths
  const double a0n = length(a0);
  const double a1n = length(a1);
  const double a2n = length(a2);

  // half cell vectors
  const double a0h = 0.5 * a0n;
  const double a1h = 0.5 * a1n;
  const double a2h = 0.5 * a2n;

  for ( int i = 0; i < nst; i++ )
    correction_[i] = D3vector(0,0,0);

  int niter = 1;
  // check if override from the debug variable
  // use: set debug mlwf_ref_niter <value>
  const map<string,string>& debug_map = s_.ctrl.debug;

  map<string,string>::const_iterator imap =
    debug_map.find("MLWF_REF_NITER");
  if ( imap != debug_map.end() )
  {
    const string val = imap->second;
    istringstream is(val);
    is >> niter;
    if ( MPIdata::onpe0() )
      cout << " override mlwf_ref_niter value: niter = " << niter << endl;
    assert(niter > 0);
  }

  for ( int iter = 0; iter < niter; iter++ )
  {
    fill(ref.begin(),ref.end(),0.0);

    // wavefunctions in real space
    vector<complex<double> > wftmp(np012loc);
    vector<complex<double> > wftmp0(np012loc);
    vector<complex<double> > wftmp1(np012loc);
    vector<complex<double> > wftmp2(np012loc);

    for ( int in = 0; in < nloc; in++ )
    {
      int n = c.jglobal(in);
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
      // position of the corrected center
      const double xc = mlwfc_[n][0] + correction_[n][0];
      const double yc = mlwfc_[n][1] + correction_[n][1];
      const double zc = mlwfc_[n][2] + correction_[n][2];

      // position of the corrected center in the a basis
      const double *bmat = cell.bmat();
      // c = B^T * (xc,yc,zc)/(2*pi)
      const double itwopi = 1.0 / ( 2.0 * M_PI );
      const double c0 = itwopi * ( bmat[0] * xc + bmat[1] * yc + bmat[2] * zc );
      const double c1 = itwopi * ( bmat[3] * xc + bmat[4] * yc + bmat[5] * zc );
      const double c2 = itwopi * ( bmat[6] * xc + bmat[7] * yc + bmat[8] * zc );

#if 0
      // compute shifted sawtooth potentials v0, v1, v2
      vector<double> v0(ft.np0()), v1(ft.np1()), v2(ft.np2());
      for ( int i = 0; i < v0.size(); i++ )
      {
        double x = i * d0 - x0;
        if ( x >  a0h ) x -= a0n;
        if ( x < -a0h ) x += a0n;
#ifdef NO_VSST
        v0[i] = x / a0n;
#else
        v0[i] = vsst(x/a0n);
#endif
      }
      for ( int j = 0; j < v1.size(); j++ )
      {
        double y = j * d1 - x1;
        if ( y >  a1h ) y -= a1n;
        if ( y < -a1h ) y += a1n;
#ifdef NO_VSST
        v1[j] = y / a1n;
#else
        v1[j] = vsst(y/a1n);
#endif
      }
      for ( int k = 0; k < v2.size(); k++ )
      {
        double z = k * d2 - x2;
        if ( z >  a2h ) z -= a2n;
        if ( z < -a2h ) z += a2n;
#ifdef NO_VSST
        v2[k] = z / a2n;
#else
        v2[k] = vsst(z/a2n);
#endif
      }
#endif

      for ( int i = 0; i < np012loc; i++ )
      {
        int i0 = ft.i(i);
        int i1 = ft.j(i);
        int i2 = ft.k(i);

        // position r in a coordinates:
        // r = i0/np0 * a0 + i1/np1 * a1 + i2/np2 * a2
        // subtract correction rc = c0 * a0 + c1 * a1 + c2 * a2
        // project (r-rc) onto reciprocal lattice vectors b_k
        // (r-rc)*b0/(2*pi) = (i0/np0-c0)
        // (r-rc)*b1/(2*pi) = (i1/np1-c1)
        // (r-rc)*b2/(2*pi) = (i2/np2-c2)

        double arg0 = static_cast<double>(i0)/ft.np0() - c0;
        double arg1 = static_cast<double>(i1)/ft.np1() - c1;
        double arg2 = static_cast<double>(i2)/ft.np2() - c2;

        // compute potential values v0, v1, v2
        // v0 = vsst((r-rc)*b0)
        // v1 = vsst((r-rc)*b1)
        // v2 = vsst((r-rc)*b2)

        // fold arguments in [-1/2,1/2] for vsst
        if ( arg0 < -0.5 ) arg0 += 1.0;
        if ( arg0 >  0.5 ) arg0 -= 1.0;
        if ( arg1 < -0.5 ) arg1 += 1.0;
        if ( arg1 >  0.5 ) arg1 -= 1.0;
        if ( arg2 < -0.5 ) arg2 += 1.0;
        if ( arg2 >  0.5 ) arg2 -= 1.0;

        const double v0 = vsst(arg0);
        const double v1 = vsst(arg1);
        const double v2 = vsst(arg2);

        // multiply wft by each potential

        const double wft = real(wftmp[i]);

        const double wft0 = v0 * wft;
        const double wft1 = v1 * wft;
        const double wft2 = v2 * wft;

        // accumulate sum for expectation values of vsst((r-rc)*b_k) k=0,1,2
        pref[0] += wft * wft0;
        pref[1] += wft * wft1;
        pref[2] += wft * wft2;

        wftmp0[i] = wft0;
        wftmp1[i] = wft1;
        wftmp2[i] = wft2;

        if ( compute_quadrupole_ )
        {
          //!! quadrupole must be computed in absolute coordinates
          throw runtime_error("quadrupole not implemented");
          pref[3] += wft0 * wft0;
          pref[4] += wft1 * wft1;
          pref[5] += wft2 * wft2;
          pref[6] += wft0 * wft1;
          pref[7] += wft1 * wft2;
          pref[8] += wft2 * wft0;
        }
      } // for i
      tmap["real"].stop();

      // ft to get xwf in reciprocal space at the last iteration
      if ( iter == niter - 1 )
      {
        tmap["ft"].start();
          ft.forward(&wftmp0[0],cx.valptr(mloc*in));
          ft.forward(&wftmp1[0],cy.valptr(mloc*in));
          ft.forward(&wftmp2[0],cz.valptr(mloc*in));
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
      for ( int ist = 0; ist < nst; ist++ )
      {
        D3vector& pcor = correction_[ist];
        D3tensor& pquad = quad_[ist];

        //!! next lines not corrected for tranformation to absolute coordinates
        throw runtime_error("ElectricEnthalpy: not implemented");
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
      // corrections in b coordinates in vector (d0,d1,d2)
      // compute correction in absolute coordinates using
      // (dx,dy,dz) = A * (d0,d1,d2)

      for ( int ist = 0; ist < nst; ist++ )
      {
        const double d0 = ref[ist*3+0]/np012v;
        const double d1 = ref[ist*3+1]/np012v;
        const double d2 = ref[ist*3+2]/np012v;
        const double dx = d0 * a0.x + d1 * a1.x + d2 * a2.x;
        const double dy = d0 * a0.y + d1 * a1.y + d2 * a2.y;
        const double dz = d0 * a0.z + d1 * a1.z + d2 * a2.z;

        D3vector& pcor = correction_[ist];
        pcor[0] += dx;
        pcor[1] += dy;
        pcor[2] += dz;
      }
    }
    tmap["real"].stop();
  } // for iter
}

////////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::print(ostream& os) const
{
  if ( pol_type_ == off ) return;

  os << fixed << right << setprecision(8);
  // print MLWF centers if pol_type_ == MLWF or MLWF_REF or MLWF_REF_Q
  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
  {
    os << "<mlwfs>" << endl;
    int nst = sd_.nst();
    os << " <mlwfset spin=\"0\" size=\"" << nst << "\">" << endl;
    for ( int i = 0; i < nst; i++ )
    {
      if ( pol_type_ == mlwf )
      {
        os << " <mlwf> <center> " << setprecision(8)
           << setw(14) << mlwfc_[i].x << " "
           << setw(14) << mlwfc_[i].y << " "
           << setw(14) << mlwfc_[i].z
           << " </center> </mlwf>"
           << endl;
      }
      else
      {
        os << " <mlwf> <center> " << setprecision(8)
           << setw(14) << mlwfc_[i].x + correction_[i].x << " "
           << setw(14) << mlwfc_[i].y + correction_[i].y << " "
           << setw(14) << mlwfc_[i].z + correction_[i].z
           << " </center> </mlwf>"
           << endl;

        if ( compute_quadrupole_ )
          os << " <quad>"
             << setw(14) << quad_[i][0] << " "
             << setw(14) << quad_[i][4] << " "
             << setw(14) << quad_[i][8] << " "
             << setw(14) << quad_[i][1] << " "
             << setw(14) << quad_[i][2] << " "
             << setw(14) << quad_[i][5]
             << " </quad>" << endl;
      }
    }
    os << " </mlwfset>" << endl;
    os << "</mlwfs>" << endl;
  }

  // print dipole
  os << setprecision(8) << fixed << right;
  os << "<dipole>\n";
  os << " <dipole_ion>   "
     << setw(14) << dipole_ion_.x << " "
     << setw(14) << dipole_ion_.y << " "
     << setw(14) << dipole_ion_.z << " </dipole_ion>\n";
  os << " <dipole_el>    "
     << setw(14) << dipole_el_.x << " "
     << setw(14) << dipole_el_.y << " "
     << setw(14) << dipole_el_.z << " </dipole_el>\n";
  os << " <dipole_total> "
     << setw(14) << dipole_total_.x << " "
     << setw(14) << dipole_total_.y << " "
     << setw(14) << dipole_total_.z << " </dipole_total>\n";
  os << "</dipole>\n";

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
    D3tensor q_total = q_ion + q_mlwf;
    //traceless quadrupole
    D3tensor q_traceless = q_total;
    q_traceless.traceless();

    os << "<quadrupole> " << endl;
    os << " <quadrupole_ion> " << endl
       << q_ion
       << " </quadrupole_ion>" << endl;
    os << " <quadrupole_el> " << endl
       << q_mlwf
       << " </quadrupole_el>" << endl;
    os << " <quadrupole_total> " << endl
       << q_total
       << " </quadrupole_total>" << endl;
    os << " <traceless_quadrupole> " << endl
       << q_traceless
       << " </traceless_quadrupole>" << endl;
    char uplo = 'u';
    D3vector eigval;
    D3tensor eigvec;
    q_traceless.syev(uplo, eigval, eigvec);
    os << " <traceless_quadrupole_eigval> " << endl
       << "  " << eigval << endl
       << " </traceless_quadrupole_eigval>" << endl;
    os << " <traceless_quadrupole_eigvec> " << endl
       << eigvec
       << " </traceless_quadrupole_eigvec>" << endl;
    os << "</quadrupole> " << endl;
  }

}

////////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::set_e_field(D3vector e_field_val)
{
  e_field_ = e_field_val;
  finite_field_ = norm2(e_field_) != 0.0;
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const ElectricEnthalpy& e )
{
  e.print(os);
  return os;
}
