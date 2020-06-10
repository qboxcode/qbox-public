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
// Species.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "Species.h"
#include "spline.h"
#include "sinft.h"
#include <cmath>
#include <cassert>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
static double simpsn ( int n, double *t )
{
  const double c0 =  17.0/48.0, c1 = 59.0/48.0,
               c2 =  43.0/48.0, c3 = 49.0/48.0;
  double sum = c0 * ( t[0] + t[n-1] ) +
               c1 * ( t[1] + t[n-2] ) +
               c2 * ( t[2] + t[n-3] ) +
               c3 * ( t[3] + t[n-4] );

  for ( int i = 4; i < n-4; i++ )
  {
    sum += t[i];
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
Species::Species(string name) : name_(name),
zval_(-1), mass_(0.0), lmax_(-1), deltar_(0.0), atomic_number_(0),
llocal_(-1), nquad_(-1), rquad_(0.0),
rcps_(0.0), uri_(""), description_("undefined"), symbol_("Undef")
{}

////////////////////////////////////////////////////////////////////////////////
bool Species::initialize(double rcpsval)
{
  rcps_ = rcpsval;
  // select appropriate initialization method
  if ( type_ == NCPP )
    return initialize_ncpp();
  else if ( type_ == SLPP )
    return initialize_slpp();
  else
    throw SpeciesInitException("potential type unknown");
}

////////////////////////////////////////////////////////////////////////////////
bool Species::initialize_ncpp()
{
  assert(description_ != "undefined");

  const double fpi = 4.0 * M_PI;

  const int np = vps_[0].size();
  if (zval_ < 0) throw SpeciesInitException("zval_ < 0");
  if (rcps_ < 0.0) throw SpeciesInitException("rcps_ < 0");
  if (mass_ < 0.0) throw SpeciesInitException("mass_ < 0");
  if (lmax_ < 0 || lmax_ > 3) throw SpeciesInitException("lmax_<0 or lmax_>3");

  if (vps_.size() < lmax_+1) throw SpeciesInitException("vps_.size < lmax_+1");

  if (llocal_ < 0 || llocal_ > lmax_)
      throw SpeciesInitException("llocal_ < 0 || llocal_ > lmax_");

  if ( nquad_ == 0 ) // KB potential
  {
    for ( int l = 0; l <= lmax_; l++ )
      if ( l != llocal_ && phi_[l].size() == 0 )
        throw SpeciesInitException("phi_[l] undefined for non-local projector");
  }

  if ( nquad_ < 0 ) throw SpeciesInitException("nquad < 0");
  if ( nquad_ > 0 && rquad_ <= 0.0 )
    throw SpeciesInitException("semilocal with rquad_ <= 0");

  // compute number of non-local projectors nlm_
  nlm_ = 0;
  nop_ = lmax_ + 1;
  for ( int l = 0; l <= lmax_; l++ )
  {
    if ( l != llocal_ )
    {
       nlm_ += 2 * l + 1;
    }
    lmap_.push_back(l);
  }

  // compute ndft_: size of radial FFT array
  // ndft_ is a power of 2 larger than ( rdftmin / deltar_ )
  // minimum outer bound in (a.u.)
  const double rdftmin = 40.0;
  assert(deltar_ > 0.0);
  ndft_ = 1;
  while ( ndft_ * deltar_ < rdftmin )
    ndft_ *= 2;

  rps_spl_.resize(ndft_);
  for ( int i = 0; i < ndft_; i++ )
    rps_spl_[i] = i * deltar_;

  vps_spl_.resize(lmax_+1);
  vps_spl2_.resize(lmax_+1);
  phi_spl_.resize(lmax_+1);
  phi_spl2_.resize(lmax_+1);

  for ( int l = 0; l <= lmax_; l++ )
  {
    vps_spl_[l].resize(ndft_);
    vps_spl2_[l].resize(ndft_);
    if ( phi_[l].size() > 0 )
    {
      phi_spl_[l].resize(ndft_);
      phi_spl2_[l].resize(ndft_);
    }
  }

  // extend rps and vps_ to full mesh (up to i==ndft_-1)

  vector<double> fint(ndft_);

  wsg_.resize(lmax_+1);
  gspl_.resize(ndft_);
  vlocg_spl_.resize(ndft_);
  vlocg_spl2_.resize(ndft_);

  vnlg_spl_.resize(lmax_+1);
  vnlg_spl2_.resize(lmax_+1);

  vector<double> vlocr(ndft_);
  vector<vector<double> > vnlr(lmax_+1);

  for ( int l = 0; l <= lmax_; l++ )
  {
    vnlr[l].resize(ndft_);
    vnlg_spl_[l].resize(ndft_+1);
    vnlg_spl2_[l].resize(ndft_+1);
  }

  // Extend vps_[l][i] up to ndft_ using -zv/r

  for ( int l = 0; l <= lmax_; l++ )
  {
    for ( int i = 0; i < np; i++ )
      vps_spl_[l][i] = vps_[l][i];
    for ( int i = np; i < ndft_; i++ )
      vps_spl_[l][i] = - zval_ / rps_spl_[i];
    if ( phi_[l].size() > 0 )
    {
      for ( int i = 0; i < np; i++ )
        phi_spl_[l][i] = phi_[l][i];
      for ( int i = np; i < ndft_; i++ )
        phi_spl_[l][i] = 0.0;
    }
  }

  // compute spline coefficients of vps_spl_ and phi_spl_
  for ( int l = 0; l <= lmax_; l++ )
  {
    const double dvdr = zval_ / (rps_spl_[ndft_-1]*rps_spl_[ndft_-1]);
    spline(ndft_,&rps_spl_[0],&vps_spl_[l][0],0.0,dvdr,0,0,&vps_spl2_[l][0]);
  }
  for ( int l = 0; l <= lmax_; l++ )
  {
    if ( l != llocal_ && phi_[l].size() > 0 )
    {
      spline(ndft_,&rps_spl_[0],&phi_spl_[l][0],0.0,0.0,0,1,&phi_spl2_[l][0]);
    }
  }

  // nonlinear core correction
  if ( has_nlcc() ) initialize_nlcc();
  ztot_ = zval_;

  // local potential: subtract the long range part due to the smeared charge
  // Next line: constant is 2/sqrt(pi)
  // math.h: # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  vlocr[0] = vps_spl_[llocal_][0] + (ztot_/rcps_) * M_2_SQRTPI;
  for ( int i = 1; i < ndft_; i++ )
  {
    vlocr[i] = vps_spl_[llocal_][i] + (ztot_/rps_spl_[i]) *
               erf( rps_spl_[i]/rcps_ );
  }

  //  Prepare the function vlocr to be used later in the Bessel transforms:
  //
  //  local potential: v(G) = 4 pi \int r^2 vloc(r) sin(Gr)/(Gr) dr
  //  -> store 4 pi r dr vps_(lmax_) in vlocr(i,is)
  //
  //  the Bessel transform is then:
  //  v(G) = 1/G \sum_r sin(Gr) vlocr

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocr[i] *= fpi * rps_spl_[i] * deltar_;
  }
  //  Local potential
  //  Compute Fourier coefficients of the local potential
  //  vlocr[i] contains 4 pi r dr vpsr(lmax_)
  //  v(G) = 4 pi \int r^2 vpsr(r) sin(Gr)/(Gr) dr
  //       = 1/G \sum_r sin(Gr) vlocr
  //
  //  v(G=0) by simpson integration
  //  v(G) = 4 pi \int r^2 vpsr(r) dr
  //       = \sum_r r vlocr
  //
  //  N.B. vlocr[i] contains 4 pi r dr (vps_(lmax_)-v_pseudocharge(r))
  //  Simpson integration up to ndft_ (V is zero beyond that point)
  //  Use fint as temporary array for integration

  for ( int i = 0; i < ndft_; i++ )
  {
    fint[i] = vlocr[i] * rps_spl_[i];
  }

  const double v0 = simpsn(ndft_,&fint[0]);
  // next line: include correction for v(g=0) to be independent of rcps
  // const double v0 = simpsn(ndft_,&fint[0]) - M_PI*zval_*rcps_*rcps_;

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocg_spl_[i] = vlocr[i];
  }

  sinft(ndft_,&vlocg_spl_[0]);
  if ( has_nlcc() ) sinft(ndft_,&nlccg_spl_[0]);

  //  Divide by g's
  gspl_[0] = 0.0;
  vlocg_spl_[0] = v0;
  if ( has_nlcc() ) nlccg_spl_[0] = zcore_;
  double fac = M_PI/(ndft_*deltar_);
  for ( int i = 1; i < ndft_; i++ )
  {
    gspl_[i] = i * fac;
    vlocg_spl_[i] /= gspl_[i];
    if ( has_nlcc() ) nlccg_spl_[i] /= gspl_[i];
  }

  //  Initialize cubic spline interpolation for local potential Vloc(G)
  //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

  spline(ndft_,&gspl_[0],&vlocg_spl_[0],0.0,0.0,0,1,&vlocg_spl2_[0]);
  if ( has_nlcc() )
    spline(ndft_,&gspl_[0],&nlccg_spl_[0],0.0,0.0,0,1,&nlccg_spl2_[0]);


  // Non-local KB projectors
  if ( nquad_ == 0 && lmax_ > 0)
  {
    for ( int l = 0; l <= lmax_; l++ )
    {
      wsg_[l] = 0.0;

      if ( l != llocal_ )
      {
        // for KB potentials, compute weights wsg[l]
        //  Compute weight wsg_[l] by integration on the linear mesh
        for ( int i = 0; i < ndft_; i++ )
        {
          double tmp = phi_spl_[l][i] * rps_spl_[i];
          fint[i] = ( vps_spl_[l][i] - vps_spl_[llocal_][i] ) * tmp * tmp;
        }
        double tmp = simpsn(ndft_,&fint[0]);
        assert(tmp != 0.0);
        // Next lines: store 1/<phi|delta_v| phi > in wsg[is][l]
        wsg_[l] = 1.0 / ( deltar_ * tmp );

        //   compute non-local projectors:
        //   w(G) = Ylm(G) i^l 4 pi \int r^2 phi_l(r) j_l(Gr) v_l(r) dr
        //   -> store 4 pi v_l(r) phi_l(r) dr in vnlr[l][i]
        //   the Bessel transform is then done by
        //   l=0: j_0(Gr) = sin(Gr)/(Gr)
        //        w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
        //   l=1: j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
        //        w_1(G) = Ylm(G) i^-1 1/G^2 \sum_r sin(Gr) vnlr -
        //                 Ylm(G) i^-1 1/G   \sum_r cos(Gr) r vnlr
        //   l=2: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
        //        w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
        //                                1/G   \sum_r sin(Gr)*r  vnlr -
        //                                3/G^2 \sum_r cos(Gr)    vnlr )

        for ( int i = 0; i < ndft_; i++ )
        {
          vnlr[l][i] = fpi * deltar_ *
            ( vps_spl_[l][i] - vps_spl_[llocal_][i] ) * phi_spl_[l][i];
        }

      }
    }

    //  compute radial Fourier transforms of vnlr

    //  Next line: vnlg_spl_ dimensioned ndft_+1 since it is passed to cosft1

    // Non-local potentials
    //
    // w(G) = Ylm(G) i^l 4 pi \int r^2 phi_l(r) j_l(Gr) v_l(r) dr
    // -> store 4 pi v_l(r) phi_l(r) dr in vnlr[l][i]
    // the Bessel transform is then done by
    // s: j_0(Gr) = sin(Gr)/(Gr)
    //    w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
    // p: j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
    //    w_1(G) = Ylm(G) i^-1 1/G^2 \sum_r sin(Gr) vnlr -
    //             Ylm(G) i^-1 1/G   \sum_r cos(Gr) r vnlr
    // d: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
    //    w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
    //                            1/G   \sum_r sin(Gr)*r  vnlr -
    //                            3/G^2 \sum_r cos(Gr)    vnlr )
    // vnlr[l][i] contains 4 pi dr phi_l(r) v_l(r)

    for ( int l = 0; l <= lmax_; l++ )
    {
      if ( l != llocal_ )
      {
        if ( l == 0 )
        {

          // s projector
          // w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
          //
          // G=0: Simpson integration up to ndft_
          // Use fint as temporary array for integration

          for ( int i = 0; i < ndft_; i++ )
          {
            fint[i] = vnlr[l][i] * rps_spl_[i] * rps_spl_[i];
          }
          const double v0 = simpsn(ndft_, &fint[0]);

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i] * rps_spl_[i];
          }

          sinft(ndft_,&vnlg_spl_[l][0]);

          vnlg_spl_[l][0] = v0;
          // Divide by g
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] /= gspl_[i];
          }

          //  Initialize cubic spline interpolation
          //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_spl_[l][0],0.0,0.0,0,1,
                 &vnlg_spl2_[l][0]);

        }
        else if ( l == 1 )
        {
          //  p projectors
          //  w_1(G) = Ylm(G) i 1/G^2 \sum_r sin(Gr) vnlr -
          //           Ylm(G) i 1/G   \sum_r cos(Gr) r vnlr
          //  vnlr(i,is,l) contains 4 pi dr phi_l(r) v_l(r)

          //  Next line: v(G=0) is zero (j_1(Gr) -> 0 as G->0 )
          const double v0 = 0.0;

          //  First part: 1/G^2 \sum_r sin(Gr) vnlr
          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i];
          }

          sinft(ndft_,&vnlg_spl_[l][0]);

          // Divide by g**2 and store in fint */
          fint[0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] = vnlg_spl_[l][i] / ( gspl_[i] * gspl_[i] );
          }

          //  Second part: cosine transform: 1/G   \sum_r cos(Gr) r vnlr

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i] * rps_spl_[i];
          }

          //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
          //  since it is used and modified by cosft1
          //  vnlg_ was dimensioned ndft_[is]+1

          vnlg_spl_[l][ndft_] = 0.0;
          cosft1(ndft_,&vnlg_spl_[l][0]);

          // Divide by g and add to fint to get vnlg_
          vnlg_spl_[l][0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = fint[i] - vnlg_spl_[l][i]/gspl_[i];
          }

          // Initialize spline interpolation
          // Use natural bc (y"=0) at G=0 and natural bc (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_spl_[l][0],0.0,0.0,1,1,
                 &vnlg_spl2_[l][0]);
        }
        else if ( l == 2 )
        {
          // d projectors
          // d: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
          //    w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
          //                            1/G   \sum_r sin(Gr)*r  vnlr -
          //                            3/G^2 \sum_r cos(Gr)    vnlr )
          // vnlr[l][i] contains 4 pi dr phi_l(r) v_l(r)

          //  Next line: v(G=0) is zero (j_2(Gr) -> 0 as G->0 )
          const double v0 = 0.0;

          // First part: sine transform 3/G^3 \sum_r sin(Gr)/r vnlr
          // Note: the integrand is linear near r=0 since vnlr(r) ~ r^2
          vnlg_spl_[l][0] = 0.0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i] / rps_spl_[i];
          }

          sinft(ndft_,&vnlg_spl_[l][0]);

          // multiply by 3/G^3 and store in fint */
          fint[0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] = 3.0 * vnlg_spl_[l][i] / ( gspl_[i]*gspl_[i]*gspl_[i] );
          }

          // Second part: sine transform -1/G \sum_r sin(Gr)*r vnlr
          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i] * rps_spl_[i];
          }

          sinft(ndft_,&vnlg_spl_[l][0]);

          // multiply by -1/G and accumulate in fint */
          fint[0] += v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] += - vnlg_spl_[l][i] / gspl_[i];
          }

          // Third part: cosine transform: -3/G^2 \sum_r cos(Gr) vnlr

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = vnlr[l][i];
          }

          //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
          //  since it is used and modified by cosft1
          //  vnlg_ was dimensioned ndft_[is]+1

          vnlg_spl_[l][ndft_] = 0.0;
          cosft1(ndft_,&vnlg_spl_[l][0]);

          // Multiply by -3/G^2 and add to fint
          fint[0] += v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] += - 3.0 * vnlg_spl_[l][i] / (gspl_[i] * gspl_[i]);
          }

          vnlg_spl_[l][0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_spl_[l][i] = fint[i];
          }

          // Initialize spline interpolation
          // Use zero first derivative at G=0 and natural (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_spl_[l][0],0.0,0.0,0,1,
                 &vnlg_spl2_[l][0]);
        }
      } // l != llocal_
    } // l
  } // nquad_ == 0 && lmax_ > 0 (KB projectors)
  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool Species::initialize_slpp()
{
  assert(description_ != "undefined");

  const double fpi = 4.0 * M_PI;

  const int np = vlocal_.size();
  nquad_ = 0;

  // number of projectors
  nop_ = 0;
  for ( int l = 0; l < proj_.size(); l++ )
    for ( int n = 0; n < proj_[l].size(); n++ )
    {
      lmap_.push_back(l);
      nop_++;
    }

  // check sanity of input
  if (zval_ < 0) throw SpeciesInitException("zval_ < 0");
  if (rcps_ < 0.0) throw SpeciesInitException("rcps_ < 0");
  if (mass_ < 0.0) throw SpeciesInitException("mass_ < 0");

  // check consistency of input
  // a) for each angular momentum l there has to be a D matrix of
  //    the size n x n, where n are the number of projectors with
  //    this angular momentum l
  if ( proj_.size() != d_.size() )
    throw SpeciesInitException("projector and dmatrix inconsistent");
  for ( int l = 0; l < proj_.size(); l++ )
  {
    // trivial case only one projector (n = 1)
    // if D matrix is not present, create a trivial one D = 1
    // otherwise leave D matrix unchanged
    if ( proj_[l].size() == 1 )
    {
      if ( d_[l].size() > 1 )
        throw SpeciesInitException("dmatrix has wrong dimension");
      d_[l].resize(1);
      if ( d_[l][0].size() > 1)
        throw SpeciesInitException("dmatrix has wrong dimension");
      d_[l][0].resize(1,1);
    }
    // general case dim( d[l] ) = n
    else
    {
      const int nmax = proj_[l].size();
      if ( d_[l].size() != nmax )
        throw SpeciesInitException("dmatrix has wrong dimension");
      for ( int n = 0; n < nmax; n++ )
      {
        if ( d_[l][n].size() != nmax )
          throw SpeciesInitException("dmatrix has wrong dimension");
      }
    }
  }

  // compute number of non-local projectors nlm_
  nlm_ = 0;
  for ( int l = 0; l < proj_.size(); l++ )
  {
    nlm_ += (2 * l + 1) * proj_[l].size();
  }

  // compute ndft_: size of radial FFT array
  // ndft_ is a power of 2 larger than ( rdftmin / deltar_ )
  // minimum outer bound in (a.u.)
  const double rdftmin = 40.0;
  // next line: limit small mesh sizes
  assert(deltar_ > 0.0001);
  ndft_ = 1;
  while ( ndft_ * deltar_ < rdftmin )
    ndft_ *= 2;

  rps_spl_.resize(ndft_);
  for ( int i = 0; i < ndft_; i++ )
    rps_spl_[i] = i * deltar_;

  phi_spl_.resize(nop_);
  phi_spl2_.resize(nop_);

  for ( int iop = 0; iop < nop_; iop++ )
  {
    phi_spl_[iop].resize(ndft_);
    phi_spl2_[iop].resize(ndft_);
  }


  // extend rps and vps_ to full mesh (up to i==ndft_-1)

  vector<double> fint(ndft_);

  gspl_.resize(ndft_);
  vlocg_spl_.resize(ndft_);
  vlocg_spl2_.resize(ndft_);

  vnlg_spl_.resize(nop_);
  vnlg_spl2_.resize(nop_);

  vector<double> vlocr(ndft_);
  vector<vector<double> > vnlr(nop_);

  for ( int iop = 0; iop < nop_; iop++ )
  {
    vnlr[iop].resize(ndft_);
    vnlg_spl_[iop].resize(ndft_ + 1);
    vnlg_spl2_[iop].resize(ndft_ + 1);
  }

  // Extend vps_[l][i] up to ndft_ using -zv/r

  // local potential
  for ( int i = 0; i < np; i++ )
    vlocr[i] = vlocal_[i];
  for ( int i = np; i < ndft_; i++ )
    vlocr[i] = -zval_ / rps_spl_[i];

  for ( int iop = 0; iop < nop_; iop++ )
  {
    // counter for projectors of same l
    int n = 0;
    if ( iop == 0 || lmap_[iop] != lmap_[iop - 1] ) n = 0;
    else n++;
    for ( int i = 0; i < np; i++ )
      phi_spl_[iop][i] = proj_[lmap_[iop]][n][i];
    for ( int i = np; i < ndft_; i++ )
      phi_spl_[iop][i] = 0.0;
  }

  for ( int iop = 0; iop < nop_; iop++ )
  {
    spline(ndft_,&rps_spl_[0],&phi_spl_[iop][0],0.0,0.0,0,1,&phi_spl2_[iop][0]);
  }

  // nonlinear core correction
  if ( has_nlcc() ) initialize_nlcc();
  ztot_ = zval_;

  // local potential: subtract the long range part due to the smeared charge
  // Next line: constant is 2/sqrt(pi)
  // math.h: # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  vlocr[0] += ( ztot_ / rcps_ ) * M_2_SQRTPI;
  for ( int i = 1; i < ndft_; i++ )
  {
    vlocr[i] += ( ztot_ / rps_spl_[i] ) * erf(rps_spl_[i] / rcps_);
  }

  //  Prepare the function vlocr to be used later in the Bessel transforms:
  //
  //  local potential: v(G) = 4 pi \int r^2 vloc(r) sin(Gr)/(Gr) dr
  //  -> store 4 pi r dr vps_(lmax_) in vlocr(i,is)
  //
  //  the Bessel transform is then:
  //  v(G) = 1/G \sum_r sin(Gr) vlocr

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocr[i] *= fpi * rps_spl_[i] * deltar_;
  }
  //  Local potential
  //  Compute Fourier coefficients of the local potential
  //  vlocr[i] contains 4 pi r dr vpsr(lmax_)
  //  v(G) = 4 pi \int r^2 vpsr(r) sin(Gr)/(Gr) dr
  //       = 1/G \sum_r sin(Gr) vlocr
  //
  //  v(G=0) by simpson integration
  //  v(G) = 4 pi \int r^2 vpsr(r) dr
  //       = \sum_r r vlocr
  //
  //  N.B. vlocr[i] contains 4 pi r dr (vps_(lmax_)-v_pseudocharge(r))
  //  Simpson integration up to ndft_ (V is zero beyond that point)
  //  Use fint as temporary array for integration

  for ( int i = 0; i < ndft_; i++ )
  {
    fint[i] = vlocr[i] * rps_spl_[i];
  }

  const double v0 = simpsn(ndft_,&fint[0]);
  // next line: include correction for v(g=0) to be independent of rcps
  // const double v0 = simpsn(ndft_,&fint[0]) - M_PI*zval_*rcps_*rcps_;

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocg_spl_[i] = vlocr[i];
  }

  sinft(ndft_,&vlocg_spl_[0]);
  if ( has_nlcc() ) sinft(ndft_,&nlccg_spl_[0]);

  //  Divide by g's
  gspl_[0] = 0.0;
  vlocg_spl_[0] = v0;
  if ( has_nlcc() ) nlccg_spl_[0] = zcore_;
  double fac = M_PI/(ndft_*deltar_);
  for ( int i = 1; i < ndft_; i++ )
  {
    gspl_[i] = i * fac;
    vlocg_spl_[i] /= gspl_[i];
    if ( has_nlcc() ) nlccg_spl_[i] /= gspl_[i];
  }

  //  Initialize cubic spline interpolation for local potential Vloc(G)
  //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

  spline(ndft_,&gspl_[0],&vlocg_spl_[0],0.0,0.0,0,1,&vlocg_spl2_[0]);
  if ( has_nlcc() )
    spline(ndft_,&gspl_[0],&nlccg_spl_[0],0.0,0.0,0,1,&nlccg_spl2_[0]);


  // Non-local projectors
  //                                          inf
  //   /                                        /
  //  |  3   i G r                     l       |     2
  //  | d r e      f (r) Y (r) = 4 pi i  Y (G) | dr r  j (Gr) f (r)
  //  |             l     lm              lm   |        l      l
  // /                                        /
  //                                          0
  // where f_l(r) = phi_l(r)

  // multiply projectors with weight
  // w_l = 1 / < phi_l | V_l | phi_l >
  // as there is no V_l all weights are 1
  wsg_.resize(nop_,1.0);
  for ( int iop = 0; iop < nop_; iop++ )
  {
    //   compute non-local projectors:
    //   w(G) = Ylm(G) i^l 4 pi \int r^2 j_l(Gr) f_l(r) dr
    //   -> store 4 pi f_l(r) dr in vnlr[l][i]
    //   the Bessel transform is then done by
    //   l=0: j_0(Gr) = sin(Gr)/(Gr)
    //        w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
    //   l=1: j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
    //        w_1(G) = Ylm(G) i^-1 1/G^2 \sum_r sin(Gr) vnlr -
    //                 Ylm(G) i^-1 1/G   \sum_r cos(Gr) r vnlr
    //   l=2: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
    //        w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
    //                                1/G   \sum_r sin(Gr)*r  vnlr -
    //                                3/G^2 \sum_r cos(Gr)    vnlr )

    for ( int i = 0; i < ndft_; i++ )
    {
      vnlr[iop][i] = fpi * deltar_ * phi_spl_[iop][i];
    }

  } // loop over projector

  //  compute radial Fourier transforms of vnlr

  //  Next line: vnlg_spl_ dimensioned ndft_+1 since it is passed to cosft1

  // Non-local potentials
  //
  // w(G) = Ylm(G) i^l 4 pi \int r^2 j_l(Gr) f_l(r) dr
  // -> store 4 pi f_l(r) dr in vnlr[l][i]
  // the Bessel transform is then done by
  // s: j_0(Gr) = sin(Gr)/(Gr)
  //    w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
  // p: j_1(Gr) = sin(Gr)/(Gr)^2 - cos(Gr)/Gr
  //    w_1(G) = Ylm(G) i^-1 1/G^2 \sum_r sin(Gr) vnlr -
  //             Ylm(G) i^-1 1/G   \sum_r cos(Gr) r vnlr
  // d: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
  //    w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
  //                            1/G   \sum_r sin(Gr)*r  vnlr -
  //                            3/G^2 \sum_r cos(Gr)    vnlr )
  // f: j_3(Gr) = sin(Gr)*(15/(Gr)^4-6/(Gr)^2)) +
  //              cos(Gr)*(1/(Gr)-15/(Gr)^3)
  //    w_3(G) = Ylm(G) i^-3 (  15/G^4 \sum_r sin(Gr)/r^2  vnlr
  //                          - 6/G^2  \sum_r sin(Gr)      vnlr
  //                          + 1/G    \sum_r cos(Gr)*r    vnlr
  //                          - 15/G^3 \sum_r cos(Gr)/r    vnlr )

  for ( int iop = 0; iop < nop_; iop++ )
  {

    if ( lmap_[iop] == 0 )
    {
      // s projector
      // w_0(G) = Ylm(G)  1/G \sum_r sin(Gr) r vnlr
      //
      // G=0: Simpson integration up to ndft_
      // Use fint as temporary array for integration

      for ( int i = 0; i < ndft_; i++ )
      {
        const double r = rps_spl_[i];
        fint[i] = vnlr[iop][i] * r * r;
      }
      const double v0 = simpsn(ndft_,&fint[0]);

      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] * rps_spl_[i];
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      vnlg_spl_[iop][0] = v0;
      // Divide by g
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] /= gspl_[i];
      }

      //  Initialize cubic spline interpolation
      //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

      spline(ndft_,&gspl_[0],&vnlg_spl_[iop][0],
             0.0,0.0,0,1,&vnlg_spl2_[iop][0]);

    }
    else if ( lmap_[iop] == 1 )
    {
      //  p projectors
      //  w_1(G) = Ylm(G) i 1/G^2 \sum_r sin(Gr) vnlr -
      //           Ylm(G) i 1/G   \sum_r cos(Gr) r vnlr

      //  Next line: v(G=0) is zero (j_1(Gr) -> 0 as G->0 )
      const double v0 = 0.0;

      //  First part: 1/G^2 \sum_r sin(Gr) vnlr
      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i];
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      // Divide by g**2 and store in fint */
      fint[0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] = vnlg_spl_[iop][i] / ( g * g );
      }

      //  Second part: cosine transform: 1/G   \sum_r cos(Gr) r vnlr

      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] * rps_spl_[i];
      }

      //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
      //  since it is used and modified by cosft1
      //  vnlg_ was dimensioned ndft_[is]+1

      vnlg_spl_[iop][ndft_] = 0.0;
      cosft1(ndft_,&vnlg_spl_[iop][0]);

      // Divide by g and add to fint to get vnlg_
      vnlg_spl_[iop][0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = fint[i] - vnlg_spl_[iop][i] / gspl_[i];
      }

      // Initialize spline interpolation
      // Use natural bc (y"=0) at G=0 and natural bc (y"=0) at Gmax

      spline(ndft_,&gspl_[0],&vnlg_spl_[iop][0],
             0.0,0.0,1,1,&vnlg_spl2_[iop][0]);
    }
    else if ( lmap_[iop] == 2 )
    {
      // d projectors
      // d: j_2(Gr) = sin(Gr)*(3/(Gr)^3 -1/(Gr)) - 3*cos(Gr)/(Gr)^2
      //    w_2(G) = Ylm(G) i^-2 (  3/G^3 \sum_r sin(Gr)/r  vnlr -
      //                            1/G   \sum_r sin(Gr)*r  vnlr -
      //                            3/G^2 \sum_r cos(Gr)    vnlr )

      //  Next line: v(G=0) is zero (j_2(Gr) -> 0 as G->0 )
      const double v0 = 0.0;

      // First part: sine transform 3/G^3 \sum_r sin(Gr)/r vnlr
      // Note: the integrand is linear near r=0 since vnlr(r) ~ r^2
      vnlg_spl_[iop][0] = 0.0;
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] / rps_spl_[i];
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      // multiply by 3/G^3 and store in fint */
      fint[0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] = 3.0 * vnlg_spl_[iop][i] / ( g * g * g );
      }

      // Second part: sine transform -1/G \sum_r sin(Gr)*r vnlr
      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] * rps_spl_[i];
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      // multiply by -1/G and accumulate in fint */
      for ( int i = 1; i < ndft_; i++ )
      {
        fint[i] += -vnlg_spl_[iop][i] / gspl_[i];
      }

      // Third part: cosine transform: -3/G^2 \sum_r cos(Gr) vnlr
      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i];
      }

      //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
      //  since it is used and modified by cosft1
      //  vnlg_ was dimensioned ndft_[is]+1

      vnlg_spl_[iop][ndft_] = 0.0;
      cosft1(ndft_,&vnlg_spl_[iop][0]);

      // Multiply by -3/G^2 and add to fint
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] += -3.0 * vnlg_spl_[iop][i] / ( g * g );
      }

      vnlg_spl_[iop][0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = fint[i];
      }

      // Initialize spline interpolation
      // Use zero first derivative at G=0 and natural (y"=0) at Gmax

      spline(ndft_,&gspl_[0],&vnlg_spl_[iop][0],
             0.0,0.0,0,1,&vnlg_spl2_[iop][0]);
    }
    else if ( lmap_[iop] == 3 )
    {
      // f projectors
      // f: j_3(Gr) = sin(Gr)*(15/(Gr)^4-6/(Gr)^2)) +
      //              cos(Gr)*(1/(Gr)-15/(Gr)^3)
      //    w_3(G) = Ylm(G) i^-3 (  15/G^4 \sum_r sin(Gr)/r^2  vnlr
      //                          - 6/G^2  \sum_r sin(Gr)      vnlr
      //                          + 1/G    \sum_r cos(Gr)*r    vnlr
      //                          - 15/G^3 \sum_r cos(Gr)/r    vnlr )

      //  Next line: v(G=0) is zero (j_3(Gr) -> 0 as G->0 )
      const double v0 = 0.0;

      // First part: sine transform 15/G^4 \sum_r sin(Gr)/r^2 vnlr
      // Note: the integrand is linear near r=0 since vnlr(r) ~ r^3
      vnlg_spl_[iop][0] = 0.0;
      for ( int i = 1; i < ndft_; i++ )
      {
        const double r = rps_spl_[i];
        vnlg_spl_[iop][i] = vnlr[iop][i] / (r*r);
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      // multiply by 15/G^4 and store in fint */
      fint[0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] = 15.0 * vnlg_spl_[iop][i] / (g*g*g*g);
      }

      // Second part: sine transform -6/G^2 \sum_r sin(Gr) vnlr
      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i];
      }

      sinft(ndft_,&vnlg_spl_[iop][0]);

      // multiply by -6/G^2 and accumulate in fint */
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] += -6.0 * vnlg_spl_[iop][i] / (g*g);
      }

      // Third part: cosine transform:  1/G \sum_r cos(Gr)*r vnlr
      for ( int i = 0; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] * rps_spl_[i];
      }

      //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
      //  since it is used and modified by cosft1
      //  vnlg_ was dimensioned ndft_[is]+1

      vnlg_spl_[iop][ndft_] = 0.0;
      cosft1(ndft_,&vnlg_spl_[iop][0]);

      // Multiply by 1/G and add to fint
      for ( int i = 1; i < ndft_; i++ )
      {
        fint[i] += vnlg_spl_[iop][i] / gspl_[i];
      }

      // Fourth part: cosine transform:  -15/G^3 \sum_r cos(Gr)/r vnlr
      vnlg_spl_[iop][0] = 0.0;
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = vnlr[iop][i] / rps_spl_[i];
      }

      //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
      //  since it is used and modified by cosft1
      //  vnlg_ was dimensioned ndft_[is]+1

      vnlg_spl_[iop][ndft_] = 0.0;
      cosft1(ndft_,&vnlg_spl_[iop][0]);

      // Multiply by -15/G^3 and add to fint
      for ( int i = 1; i < ndft_; i++ )
      {
        const double g = gspl_[i];
        fint[i] += -15.0 * vnlg_spl_[iop][i] / (g*g*g);
      }

      vnlg_spl_[iop][0] = v0;
      for ( int i = 1; i < ndft_; i++ )
      {
        vnlg_spl_[iop][i] = fint[i];
      }

      // Initialize spline interpolation
      // Use zero first derivative at G=0 and natural (y"=0) at Gmax

      spline(ndft_,&gspl_[0],&vnlg_spl_[iop][0],
             0.0,0.0,0,1,&vnlg_spl2_[iop][0]);
    }
    else
    {
      throw SpeciesInitException("Species::initialize_slpp: l < 0 or l > 3");
    }
  } // loop over projector
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void Species::initialize_nlcc()
{
  const double fpi = 4.0 * M_PI;

  // allocate array
  nlcc_spl_.resize(ndft_);
  nlcc_spl2_.resize(ndft_);
  nlccg_spl_.resize(ndft_);
  nlccg_spl2_.resize(ndft_);

  fill( nlcc_spl_.begin(), nlcc_spl_.end(), 0.0 );
  copy( nlcc_.begin(), nlcc_.end(), nlcc_spl_.begin() );

  // compute spline coefficients
  spline(ndft_,&rps_spl_[0],&nlcc_spl_[0],0.0,0.0,0,1,&nlcc_spl2_[0]);

  // prepare the function nlccr to be used later in the Bessel transforms
  //
  //  local density: nlcc(G) = 4 pi \int r^2 nlcc(r) sin(Gr)/(Gr) dr
  //  -> store 4 pi r dr nlcc_spl_ in nlccr(i)
  //
  //  the Bessel transform is then:
  //  nlcc(G) = 1/G \sum_r sin(Gr) nlccr

  for ( int i = 0; i < ndft_; i++ )
  {
    nlccg_spl_[i] = fpi * rps_spl_[i] * deltar_ * nlcc_spl_[i];
  }

  // integrate core correction density
  //   /         2
  //  | dr 4 pi r n(r)
  // /
  // nlccr already contains factor 4 pi r
  // use fint for integration
  vector<double> fint(ndft_);
  for ( int i = 0; i < ndft_; i++ )
  {
    fint[i] = nlccg_spl_[i] * rps_spl_[i];
  }
  zcore_ = simpsn(ndft_,&fint[0]);

}

void Species::phi(int l, double r, double &val)
{
  val = 0.0;
  if ( phi_[l].size() == 0 || l > lmax_ || r > rps_spl_[ndft_-1] )
    return;
  splint(ndft_,&rps_spl_[0],&phi_spl_[l][0],&phi_spl2_[l][0],r,&val);
}

void Species::vpsr(int l, double r, double &v)
{
  if ( l > lmax_ || r > rps_spl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(ndft_,&rps_spl_[0],&vps_spl_[l][0],&vps_spl2_[l][0],r,&v);
  }
}

void Species::dvpsr(int l, double r, double &v, double &dv)
{
  if ( l > lmax_ || r > rps_spl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(ndft_,&rps_spl_[0],&vps_spl_[l][0],&vps_spl2_[l][0],r,&v,&dv);
  }
}

void Species::vlocg(double g, double &v)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(ndft_,&gspl_[0],&vlocg_spl_[0],&vlocg_spl2_[0],g,&v);
  }
}

void Species::dvlocg(double g, double &v, double &dv)
{
  if ( g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(ndft_,&gspl_[0],&vlocg_spl_[0],&vlocg_spl2_[0],g,&v,&dv);
  }
}

void Species::vnlg(int iop, double g, double &v)
{
  assert ( iop >= 0 && iop < nop_ );
  int l = lmap_[iop];
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(ndft_,&gspl_[0],&vnlg_spl_[iop][0],&vnlg_spl2_[iop][0],g,&v);
  }
}

void Species::dvnlg(int iop, double g, double &v, double &dv)
{
  assert ( iop >= 0 && iop < nop_ );
  int l = lmap_[iop];
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(ndft_,&gspl_[0],&vnlg_spl_[iop][0],&vnlg_spl2_[iop][0],g,&v,&dv);
  }
}

double Species::rhopsg( double g )
{
  double arg = 0.25 * rcps_ * rcps_ * g * g;
  return -ztot_ * exp( -arg );
}

// core density for nonlinear core correction in reciprocal space
void Species::rhocoreg(double g, double &rho)
{
  if ( has_nlcc() && g <= gspl_[ndft_-1] )
  {
    // spline interpolation
    splint(ndft_,&gspl_[0],&nlccg_spl_[0],&nlccg_spl2_[0],g,&rho);
  }
  else
  {
    rho = 0.0;
  }
}

// core density for nonlinear core correction in reciprocal space
void Species::drhocoreg(double g, double &rho, double &drho)
{
  if ( has_nlcc() && g <= gspl_[ndft_ - 1] )
  {
    // spline interpolation
    splintd(ndft_,&gspl_[0],&nlccg_spl_[0],&nlccg_spl2_[0],g,&rho,&drho);
  }
  else
  {
    rho = 0.0;
    drho = 0.0;
  }
}

ostream& operator << ( ostream &os, Species &s )
{
  s.print(os, false);
  return os;
}

void Species::print(ostream &os, bool expanded_form)
{
  // XML output of species
  // print in compact form if the uri is known
  // and if expanded_form==false
  if ( !expanded_form && uri() != "" )
  {
    os <<"<species name=\"" << name()
       << "\" href=\"" << uri() << "\"/>" << endl;
  }
  else
  {
    os.setf(ios::scientific,ios::floatfield);
    os << setprecision(12);
    os <<"<species name=\"" << name() << "\">" << endl;
    os << "<description>" << description() << "</description>" << endl;
    os << "<symbol>" << symbol() << "</symbol>" << endl;
    os << "<atomic_number>" << atomic_number() << "</atomic_number>" << endl;
    os << "<mass>" << mass() << "</mass>" << endl;
    if ( type_ == NCPP )
    {
      os << "<norm_conserving_pseudopotential>" << endl;
      os << "<valence_charge>" << zval() << "</valence_charge>" << endl;
      os << "<lmax>" << lmax() << "</lmax>" << endl;
      os << "<llocal>" << llocal() << "</llocal>" << endl;
      os << "<nquad>" << nquad() << "</nquad>" << endl;
      os << "<rquad>" << rquad() << "</rquad>" << endl;
      os << "<mesh_spacing>" << deltar() << "</mesh_spacing>" << endl;
      if ( nlcc_.size() > 0 ) print_nlcc(os);
      for ( int l = 0; l <= lmax(); l++ )
      {
        const int size = vps_[l].size();
        os << "<projector l=\"" << l << "\" size=\"" << size
           << "\">" << endl;
        os << "<radial_potential>\n";
        for ( int i = 0; i < size; i++ )
          os << vps_[l][i] << "\n";
        os << "</radial_potential>\n";
        os << "<radial_function>\n";
        if ( l < phi_.size() && phi_[l].size() == size )
          print_radial_function(os,phi_[l]);
        os << "</radial_function>\n";
        os << "</projector>" << endl;
      }
      os << "</norm_conserving_pseudopotential>" << endl;
    }
    else if ( type_ == SLPP )
    {
      os << "<norm_conserving_semilocal_pseudopotential>" << endl;
      os << "<valence_charge>" << zval() << "</valence_charge>" << endl;
      os << "<mesh_spacing>" << deltar() << "</mesh_spacing>" << endl;
      if ( nlcc_.size() > 0 ) print_nlcc(os);
      os << "<local_potential size=\"" << vlocal_.size() << "\">" << endl;
      print_radial_function(os,vlocal_);
      os << "</local_potential>" << endl;
      // print projector
      for ( int l = 0; l < proj_.size(); l++ )
      {
        for ( int n = 0; n < proj_[l].size(); n++ )
        {
          os << "<projector l=\"" << l << "\" i=\"" << n + 1 << "\" size=\""
            << proj_[l][n].size() << "\">" << endl;
          // print radial function
          print_radial_function(os,proj_[l][n]);
          os << "</projector>" << endl;
        }
      }
      // print D matrix
      for ( int l = 0; l < d_.size(); l++ )
      {
        for ( int n = 0; n < d_[l].size(); n++ )
          for ( int m = 0; m < d_[l][n].size(); m++ )
            os << "<d_ij l=\"" << l << "\" i=\"" << n + 1 << "\" j=\"" << m + 1
              << "\">" << d_[l][n][m] << "</d_ij>" << endl;
      }
      os << "</norm_conserving_semilocal_pseudopotential>" << endl;
    }
    os << "</species>" << endl;
  }
}

// print core density of nonlinear core correction
void Species::print_nlcc(std::ostream& os)
{
  os << "<core_density size=\"" << nlcc_.size() << "\">" << endl;
  print_radial_function(os,nlcc_);
  os << "</core_density>" << endl;
}

// print any radial function
void Species::print_radial_function(std::ostream& os,
  const std::vector<double>& rad_func)
{
  for ( int i = 0; i < rad_func.size(); i++ )
    os << rad_func[i] << endl;
}

void Species::info(ostream &os)
{
  os.setf(ios::left,ios::adjustfield);

  if ( uri() != "" )
  {
    os <<"<species name=\"" << name()
       << "\" href=\"" << uri() << "\">" << endl;
  }
  else
  {
    os <<"<species name=\"" << name() << "\">" << endl;
  }
  os << " <description>" << description() << " </description>" << endl;
  os << " <symbol>" << symbol() << "</symbol>" << endl;
  os << " <atomic_number>" << atomic_number() << "</atomic_number>" << endl;
  os << " <mass>" << mass() << "</mass>" << endl;
  if ( type_ == NCPP )
  {
    os << " <norm_conserving_pseudopotential>" << endl;
    os << " <valence_charge>" << zval() << "</valence_charge>" << endl;
    os << " <lmax>" << lmax() << "</lmax>" << endl;
    os << " <llocal>" << llocal() << "</llocal>" << endl;
    os << " <nquad>" << nquad() << "</nquad>" << endl;
    os << " <rquad>" << rquad() << "</rquad>" << endl;
    os << " <mesh_spacing>" << deltar() << "</mesh_spacing>" << endl;
    os << " </norm_conserving_pseudopotential>" << endl;
  }
  else
  {
    os << " <norm_conserving_semilocal_pseudopotential>" << endl;
    os << " <valence_charge>" << zval() << "</valence_charge>" << endl;
    os << " <mesh_spacing>" << deltar() << "</mesh_spacing>" << endl;
    os << " </norm_conserving_semilocal_pseudopotential>" << endl;
  }
  os << "</species>" << endl;

  // describe type of potential
  if ( type_ == NCPP )
  {
    if ( nquad() == 0 )
    {
      if ( lmax() == 0 )
        os << " local potential" << endl;
      else
        os << " Kleinman-Bylander potential" << endl;
    }
    else if ( nquad() > 0 )
    {
      os << " NCPP potential with " << nquad()
         << " quadrature points in [0.0, " << rquad() << "]" << endl;
      os << " local (within 1.e-6) beyond r = " << rcut_loc(1.e-6) << endl;
      os << " local (within 1.e-5) beyond r = " << rcut_loc(1.e-5) << endl;
      os << " local (within 1.e-4) beyond r = " << rcut_loc(1.e-4) << endl;
      os << " local (within 1.e-3) beyond r = " << rcut_loc(1.e-3) << endl;
    }
  }
  else if ( type_ == SLPP )
  {
    os << " SLPP semilocal potential" << endl;
  }
  os << " rcps_ =   " << rcps() << endl;
  os.setf(ios::right,ios::adjustfield);
}

////////////////////////////////////////////////////////////////////////////////
double Species::rcut_loc(double epsilon)
{
  // find radius at which the largest deviation delta_vnl(r) is < epsilon
  double delta = 0.0;
  int i = ndft_-1;
  while ( ( delta < epsilon ) && i > 0 )
  {
    i--;
    for ( int l = 0; l <= lmax_; l++ )
    {
      if ( l != llocal_ )
      {
        // compute deviation vps_[l][i]
        double dv = fabs( vps_spl_[l][i] - vps_spl_[llocal_][i] );
        delta = dv > delta ? dv : delta;
      }
    }
  }
  // adjust i so that delta_v[i] < epsilon
  if ( i < ndft_-1 ) i++;

  return rps_spl_[i];
}

////////////////////////////////////////////////////////////////////////////////
// helper function that extract l and n from projector index
Species::ProjectorData Species::get_proj_data(int ipr)
{
  ProjectorData res;
  res.l = -1;
  res.n = 0;
  int jpr = 0;
  for ( int iop = 0; iop < nop_; iop++ )
  {
    int l = lmap_[iop];
    // n labels differentiates between projectors of same l
    if ( res.l == l ) res.n++;
    // reset n if new angular momentum is found
    else
    {
      res.l = l;
      res.n = 0;
    }
    // there is one projector Y_lm
    for ( res.m = -l; res.m <= l; res.m++ )
    {
      if ( ipr == jpr ) return res;
      jpr++;
    }
  }
  // projector out of bounds
  assert(false);
  res.l = -1;
  res.n = 0;
  return res;
}

////////////////////////////////////////////////////////////////////////////////
// extract D matrix
double Species::dmatrix(int ipr, int jpr)
{
  // extract angular momentum of projectors
  const ProjectorData iem = get_proj_data(ipr);
  const ProjectorData jem = get_proj_data(jpr);

  // projectors have only nonzero overlap if diagonal in l and m
  if ( iem.l != jem.l || iem.m != jem.m ) return 0;

  // if dmatrix is not present identity matrix is assumed
  if ( d_.size() <= iem.l || d_[iem.l].size() == 0 ) return 1;

  // return dmatrix
  return d_[iem.l][iem.n][jem.n];
}
