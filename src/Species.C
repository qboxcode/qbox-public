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
// Species.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Species.C,v 1.16 2010-02-20 23:25:38 fgygi Exp $

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
Species::Species(const Context& ctxt, string name) : ctxt_(ctxt), name_(name),
zval_(-1), mass_(0.0), lmax_(-1), deltar_(0.0), atomic_number_(0),
llocal_(-1), nquad_(-1), rquad_(0.0),
rcps_(0.0), uri_(""), description_("undefined"), symbol_("Undef")
{}

////////////////////////////////////////////////////////////////////////////////
bool Species::initialize(double rcpsval)
{
  // initialize the Species
  rcps_ = rcpsval;

  assert(description_ != "undefined");

  const double fpi = 4.0 * M_PI;

  const int np = vps_[0].size();
  if (zval_ < 0) throw SpeciesInitException("zval_ < 0");
  if (rcps_ < 0.0) throw SpeciesInitException("rcps_ < 0");
  if (mass_ < 0.0) throw SpeciesInitException("mass_ < 0");
  if (lmax_ < 0 || lmax_ > 3) throw SpeciesInitException("lmax_ <0 or lmax_ >3");

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
  for ( int l = 0; l <= lmax_; l++ )
  {
    if ( l != llocal_ )
    {
       nlm_ += 2 * l + 1;
    }
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

  rps_.resize(ndft_);
  for ( int i = 0; i < ndft_; i++ )
    rps_[i] = i * deltar_;

  vps_spl_.resize(lmax_+1);
  phi_spl_.resize(lmax_+1);

  for ( int l = 0; l <= lmax_; l++ )
  {
    vps_[l].resize(ndft_);
    phi_[l].resize(ndft_);
    vps_spl_[l].resize(ndft_);
    phi_spl_[l].resize(ndft_);
  }

  // extend rps and vps_ to full mesh (up to i==ndft_-1)

  vector<double> fint(ndft_);

  wsg_.resize(lmax_+1);
  gspl_.resize(ndft_);
  vlocg_.resize(ndft_);
  vlocg_spl.resize(ndft_);

  vnlg_.resize(lmax_+1);
  vnlg_spl.resize(lmax_+1);

  vector<double> vlocr(ndft_);
  vector<vector<double> > vnlr(lmax_+1);

  for ( int l = 0; l <= lmax_; l++ )
  {
    vnlr[l].resize(ndft_);
    vnlg_[l].resize(ndft_+1);
    vnlg_spl[l].resize(ndft_+1);
  }

  // Extend vps_[l][i] up to ndft_ using -zv/r

  for ( int l = 0; l <= lmax_; l++ )
  {
    for ( int i = np; i < ndft_; i++ )
    {
      vps_[l][i] = - zval_ / rps_[i];
    }
  }

  // compute spline coefficients of vps_ and phi_
  for ( int l = 0; l <= lmax_; l++ )
  {
    const double dvdr = zval_ / (rps_[ndft_-1]*rps_[ndft_-1]);
    spline(ndft_,&rps_[0],&vps_[l][0],0.0,dvdr,0,0,&vps_spl_[l][0]);
  }
  for ( int l = 0; l <= lmax_; l++ )
  {
    if ( l != llocal_ )
    {
      spline(ndft_,&rps_[0],&phi_[l][0],0.0,0.0,0,1,&phi_spl_[l][0]);
    }
  }

  // local potential: subtract the long range part due to the smeared charge
  // Next line: constant is 2/sqrt(pi)
  // math.h: # define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
  vlocr[0] = vps_[llocal_][0] + (zval_/rcps_) * M_2_SQRTPI;
  for ( int i = 1; i < ndft_; i++ )
  {
    vlocr[i] = vps_[llocal_][i] + (zval_/rps_[i]) * erf( rps_[i]/rcps_ );
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
    vlocr[i] *= fpi * rps_[i] * deltar_;
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
    fint[i] = vlocr[i] * rps_[i];
  }

  const double v0 = simpsn(ndft_,&fint[0]);
  // next line: include correction for v(g=0) to be independent of rcps
  // const double v0 = simpsn(ndft_,&fint[0]) - M_PI*zval_*rcps_*rcps_;

  for ( int i = 0; i < ndft_; i++ )
  {
    vlocg_[i] = vlocr[i];
  }

  sinft(ndft_,&vlocg_[0]);

  //  Divide by g's
  gspl_[0] = 0.0;
  vlocg_[0] = v0;
  double fac = M_PI/(ndft_*deltar_);
  for ( int i = 1; i < ndft_; i++ )
  {
    gspl_[i] = i * fac;
    vlocg_[i] /= gspl_[i];
  }

  //  Initialize cubic spline interpolation for local potential Vloc(G)
  //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

  spline(ndft_,&gspl_[0],&vlocg_[0],0.0,0.0,0,1,&vlocg_spl[0]);

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
          double tmp = phi_[l][i] * rps_[i];
          fint[i] = ( vps_[l][i] - vps_[llocal_][i] ) * tmp * tmp;
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
                       ( vps_[l][i] - vps_[llocal_][i] ) * phi_[l][i];
        }

      }
    }

    //  compute radial Fourier transforms of vnlr

    //  Next line: vnlg_ is dimensioned ndft_+1 since it is passed to cosft1

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
            fint[i] = vnlr[l][i] * rps_[i] * rps_[i];
          }
          const double v0 = simpsn(ndft_, &fint[0]);

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_[l][i] = vnlr[l][i] * rps_[i];
          }

          sinft(ndft_,&vnlg_[l][0]);

          vnlg_[l][0] = v0;
          // Divide by g
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_[l][i] /= gspl_[i];
          }

          //  Initialize cubic spline interpolation
          //  Use zero first derivative at G=0 and natural (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_[l][0],0.0,0.0,0,1,&vnlg_spl[l][0]);

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
            vnlg_[l][i] = vnlr[l][i];
          }

          sinft(ndft_,&vnlg_[l][0]);

          // Divide by g**2 and store in fint */
          fint[0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] = vnlg_[l][i] / ( gspl_[i] * gspl_[i] );
          }

          //  Second part: cosine transform: 1/G   \sum_r cos(Gr) r vnlr

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_[l][i] = vnlr[l][i] * rps_[i];
          }

          //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
          //  since it is used and modified by cosft1
          //  vnlg_ was dimensioned ndft_[is]+1

          vnlg_[l][ndft_] = 0.0;
          cosft1(ndft_,&vnlg_[l][0]);

          // Divide by g and add to fint to get vnlg_
          vnlg_[l][0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_[l][i] = fint[i] - vnlg_[l][i]/gspl_[i];
          }

          // Initialize spline interpolation
          // Use natural bc (y"=0) at G=0 and natural bc (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_[l][0],0.0,0.0,1,1,&vnlg_spl[l][0]);
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
          vnlg_[l][0] = 0.0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_[l][i] = vnlr[l][i] / rps_[i];
          }

          sinft(ndft_,&vnlg_[l][0]);

          // multiply by 3/G^3 and store in fint */
          fint[0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] = 3.0 * vnlg_[l][i] / ( gspl_[i] * gspl_[i] * gspl_[i] );
          }

          // Second part: sine transform -1/G \sum_r sin(Gr)*r vnlr
          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_[l][i] = vnlr[l][i] * rps_[i];
          }

          sinft(ndft_,&vnlg_[l][0]);

          // multiply by -1/G and accumulate in fint */
          fint[0] += v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] += - vnlg_[l][i] / gspl_[i];
          }

          // Third part: cosine transform: -3/G^2 \sum_r cos(Gr) vnlr

          for ( int i = 0; i < ndft_; i++ )
          {
            vnlg_[l][i] = vnlr[l][i];
          }

          //  N.B. Next line: Initialize also vnlg_[l][ndft_] to zero
          //  since it is used and modified by cosft1
          //  vnlg_ was dimensioned ndft_[is]+1

          vnlg_[l][ndft_] = 0.0;
          cosft1(ndft_,&vnlg_[l][0]);

          // Multiply by -3/G^2 and add to fint
          fint[0] += v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            fint[i] += - 3.0 * vnlg_[l][i] / (gspl_[i] * gspl_[i]);
          }

          vnlg_[l][0] = v0;
          for ( int i = 1; i < ndft_; i++ )
          {
            vnlg_[l][i] = fint[i];
          }

          // Initialize spline interpolation
          // Use zero first derivative at G=0 and natural (y"=0) at Gmax

          spline(ndft_,&gspl_[0],&vnlg_[l][0],0.0,0.0,0,1,&vnlg_spl[l][0]);
        }
      } // l != llocal_
    } // l
  } // nquad_ == 0 && lmax_ > 0 (KB projectors)
  return true;
}

void Species::phi(int l, double r, double &val)
{
  if ( l > lmax_ || r > rps_[ndft_-1] )
  {
    val = 0.0;
  }
  else
  {
    splint(ndft_,&rps_[0],&phi_[l][0],&phi_spl_[l][0],r,&val);
  }
}

void Species::vpsr(int l, double r, double &v)
{
  if ( l > lmax_ || r > rps_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(ndft_,&rps_[0],&vps_[l][0],&vps_spl_[l][0],r,&v);
  }
}

void Species::dvpsr(int l, double r, double &v, double &dv)
{
  if ( l > lmax_ || r > rps_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(ndft_,&rps_[0],&vps_[l][0],&vps_spl_[l][0],r,&v,&dv);
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
    splint(ndft_,&gspl_[0],&vlocg_[0],&vlocg_spl[0],g,&v);
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
    splintd(ndft_,&gspl_[0],&vlocg_[0],&vlocg_spl[0],g,&v,&dv);
  }
}

void Species::vnlg(int l, double g, double &v)
{
  assert ( l >= 0 && l <= lmax_ );
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
  }
  else
  {
    splint(ndft_,&gspl_[0],&vnlg_[l][0],&vnlg_spl[l][0],g,&v);
  }
}

void Species::dvnlg(int l, double g, double &v, double &dv)
{
  assert ( l >= 0 && l <= lmax_ );
  if ( l == llocal_ || g > gspl_[ndft_-1] )
  {
    v = 0.0;
    dv = 0.0;
  }
  else
  {
    splintd(ndft_,&gspl_[0],&vnlg_[l][0],&vnlg_spl[l][0],g,&v,&dv);
  }
}

double Species::rhopsg( double g )
{
  double arg = 0.25 * rcps_ * rcps_ * g * g;
  return -zval_ * exp( -arg );
}

ostream& operator << ( ostream &os, Species &s )
{
  // XML output of species
  // If the uri is known, use href to refer to it
  if ( s.uri() != "" )
  {
    os <<"<species name=\"" << s.name()
       << "\" href=\"" << s.uri() << "\"/>" << endl;
  }
  else
  {
    os <<"<species name=\"" << s.name() << "\">" << endl;
    os << "<description>" << s.description() << "</description>" << endl;
    os << "<symbol>" << s.symbol() << "</symbol>" << endl;
    os << "<atomic_number>" << s.atomic_number() << "</atomic_number>" << endl;
    os << "<mass>" << s.mass() << "</mass>" << endl;
    os << "<norm_conserving_pseudopotential>" << endl;
    os << "<valence_charge>" << s.zval() << "</valence_charge>" << endl;
    os << "<lmax>" << s.lmax() << "</lmax>" << endl;
    os << "<llocal>" << s.llocal() << "</llocal>" << endl;
    os << "<nquad>" << s.nquad() << "</nquad>" << endl;
    os << "<rquad>" << s.rquad() << "</rquad>" << endl;
    os << "<mesh_spacing>" << s.deltar() << "</mesh_spacing>" << endl;
    os.setf(ios::fixed,ios::floatfield);
    os << setprecision(6);
    for ( int l = 0; l <= s.lmax(); l++ )
    {
      const int size = s.vps()[l].size();
      os << "<projector l=\"" << l << "\" size=\"" << size
         << "\">" << endl;
      os << "<radial_potential>\n";
      for ( int i = 0; i < size; i++ )
        os << s.vps()[l][i] << "\n";
      os << "</radial_potential>\n";
      if ( l < s.phi().size() && s.phi()[l].size() == size )
      {
        os << "<radial_function>\n";
        for ( int i = 0; i < size; i++ )
          os << s.phi()[l][i] << "\n";
        os << "</radial_function>\n";
      }
      os << "</projector>" << endl;
    }
    os << "</norm_conserving_pseudopotential>" << endl;
    os << "</species>" << endl;
  }

  return os;
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
  os << " <description>" << description() << "</description>" << endl;
  os << " <symbol>" << symbol() << "</symbol>" << endl;
  os << " <atomic_number>" << atomic_number() << "</atomic_number>" << endl;
  os << " <mass>" << mass() << "</mass>" << endl;
  os << " <norm_conserving_pseudopotential>" << endl;
  os << " <valence_charge>" << zval() << "</valence_charge>" << endl;
  os << " <lmax>" << lmax() << "</lmax>" << endl;
  os << " <llocal>" << llocal() << "</llocal>" << endl;
  os << " <nquad>" << nquad() << "</nquad>" << endl;
  os << " <rquad>" << rquad() << "</rquad>" << endl;
  os << " <mesh_spacing>" << deltar() << "</mesh_spacing>" << endl;
  os << " </norm_conserving_pseudopotential>" << endl;
  os << "</species>" << endl;

  // describe type of potential
  os.setf(ios::fixed,ios::floatfield);
  os << setprecision(6);
  if ( nquad() == 0 )
  {
    if ( lmax() == 0 )
      os << " local potential" << endl;
    else
      os << " Kleinman-Bylander potential" << endl;
  }
  else
  {
    os << " Semi-local potential with " << nquad()
       << " quadrature points in [0.0, " << rquad() << "]" << endl;
    os << " local (within 1.e-6) beyond r = " << rcut_loc(1.e-6) << endl;
    os << " local (within 1.e-5) beyond r = " << rcut_loc(1.e-5) << endl;
    os << " local (within 1.e-4) beyond r = " << rcut_loc(1.e-4) << endl;
    os << " local (within 1.e-3) beyond r = " << rcut_loc(1.e-3) << endl;
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
        double dv = fabs( vps_[l][i] - vps_[llocal_][i] );
        delta = dv > delta ? dv : delta;
      }
    }
  }
  // adjust i so that delta_v[i] < epsilon
  if ( i < ndft_-1 ) i++;

  return rps_[i];
}
