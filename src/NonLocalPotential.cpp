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
// NonLocalPotential.C
//
////////////////////////////////////////////////////////////////////////////////

#include "NonLocalPotential.h"
#include "Species.h"
#include "blas.h"
#include <iomanip>
using namespace std;

#if USE_MASSV
extern "C" void vsincos(double *x, double *y, double *z, int *n);
#endif


////////////////////////////////////////////////////////////////////////////////
NonLocalPotential::~NonLocalPotential(void)
{
#if TIMING
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    if ( ctxt_.myproc()==0 )
    {
      string s = "name=\"" + (*i).first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::init(void)
{
  const int ngwl = basis_.localsize();

  nsp = atoms_.nsp();

  nop.resize(nsp);
  lloc.resize(nsp);
  lproj.resize(nsp);
  na.resize(nsp);
  npr.resize(nsp);
  nprna.resize(nsp);
  wt.resize(nsp);
  twnl.resize(nsp);
  dtwnl.resize(nsp);

  nquad.resize(nsp);
  rquad.resize(nsp);
  wquad.resize(nsp);

  nspnl = 0;
  for ( int is = 0; is < nsp; is++ )
  {
    Species *s = atoms_.species_list[is];

    npr[is] = 0;
    nprna[is] = 0;

    if ( s->non_local() )
    {
      nspnl++;
      na[is] = atoms_.na(is);
      nop[is] = s->nop();
      lloc[is] = s->llocal();
      nquad[is] = s->nquad();

      // compute total number of projectors npr[is]
      // KB potentials have nlm projectors
      // semilocal potentials have nlm*nquad projectors
      if ( s->nquad() == 0 )
      {
        npr[is] = s->nlm();
      }
      else
      {
        npr[is] = s->nlm() * nquad[is];
      }
      nprna[is] = npr[is] * na[is];

      // l value for projector ipr
      lproj[is].resize(npr[is]);

      twnl[is].resize(npr[is]*ngwl);
      dtwnl[is].resize(npr[is]*6*ngwl);

      // quadrature abcissae and weights
      rquad[is].resize(nquad[is]);
      wquad[is].resize(nquad[is]);

      enum quadrature_rule_type { TRAPEZOID, MODIF_TRAPEZOID,
      TRAPEZOID_WITH_RIGHT_ENDPOINT,
      SIMPSON };

      const quadrature_rule_type quad_rule = TRAPEZOID;
      //const quadrature_rule_type quad_rule = MODIF_TRAPEZOID;
      //const quadrature_rule_type quad_rule = TRAPEZOID_WITH_RIGHT_ENDPOINT;
      //const quadrature_rule_type quad_rule = SIMPSON;

      if ( quad_rule == TRAPEZOID )
      {
        // trapezoidal rule with interior points only
        // (end points are zero)
        const double h = s->rquad() / (nquad[is]+1);
        for ( int iquad = 0; iquad < nquad[is]; iquad++ )
        {
          rquad[is][iquad] = (iquad+1) * h;
          wquad[is][iquad] = h;
        }
        //cout << " NonLocalPotential::init: trapezoidal rule (interior)"
        //     << endl;
      }
      else if ( quad_rule == MODIF_TRAPEZOID )
      {
        // use modified trapezoidal rule with interior points, and include
        // correction for first derivative at r=0 as
        // h^2/12 f'(0) where f'(0) is estimated with f(h)/h
        // i.e. add correction h/12) * f(h)
        // See Davis & Rabinowitz, p. 132
        const double h = s->rquad() / (nquad[is]+1);
        for ( int iquad = 0; iquad < nquad[is]; iquad++ )
        {
          rquad[is][iquad] = (iquad+1) * h;
          wquad[is][iquad] = h;
        }
        wquad[is][0] += h / 12.0;
        //cout << " NonLocalPotential::init: modified trapezoidal rule"
        //     << endl;
      }
      else if ( quad_rule == TRAPEZOID_WITH_RIGHT_ENDPOINT )
      {
        const double h = s->rquad() / nquad[is];
        for ( int iquad = 0; iquad < nquad[is]; iquad++ )
        {
          rquad[is][iquad] = (iquad+1) * h;
          wquad[is][iquad] = h;
        }
        wquad[is][nquad[is]-1] = 0.5 * h;
        //cout << " NonLocalPotential::init: trapezoidal rule, right endpoint"
        //     << endl;
      }
      else if ( quad_rule == SIMPSON )
      {
        // must have 2n+1 points
        assert(nquad[is]%2==1);
        const double h = s->rquad() / (nquad[is]-1);
        for ( int iquad = 0; iquad < nquad[is]; iquad++ )
        {
          rquad[is][iquad] = iquad * h;
          if ( ( iquad == 0 ) || ( iquad == nquad[is]-1 ) )
            wquad[is][iquad] = h / 3.0;
          else if ( iquad % 2 == 0 )
            wquad[is][iquad] = h * 2.0 / 3.0;
          else
            wquad[is][iquad] = h * 4.0 / 3.0;
        }
        //cout << " NonLocalPotential::init: Simpson rule" << endl;
      }
      else
      {
        assert(false);
      }

      // compute weights wt[is][ipr]
      wt[is].resize(npr[is]);

      // compute lproj[is][ipr]
      int ipr_base = 0;
      for ( int iop = 0; iop < nop[is]; iop++ )
      {
        int l = s->l(iop);
        if ( l != lloc[is] )
        {
          if ( nquad[is] == 0 )
          {
            // Kleinman-Bylander form
            // wt[is][ipr]
            // index = ipr_base+m
            for ( int m = 0; m < 2*l+1; m++ )
            {
              const int ipr = ipr_base + m;
              wt[is][ipr] = s->wsg(iop);
              lproj[is][ipr] = l;
            }
            ipr_base += 2*l+1;
          }
          else
          {
            for ( int iquad = 0; iquad < nquad[is]; iquad++ )
            {
              const double r = rquad[is][iquad];
              double v,dv,vl,dvl;
              s->dvpsr(l,r,v,dv);
              s->dvpsr(lloc[is],r,vl,dvl);
              // wt[is][iquad+ipr*nquad]
              for ( int m = 0; m < 2*l+1; m++ )
              {
                const int ipr = ipr_base + iquad + nquad[is] * m;
                wt[is][ipr] = ( v - vl ) * wquad[is][iquad];
                lproj[is][ipr] = l;
              }
            }
            ipr_base += (2*l+1) * nquad[is];
          }
        }
      }
      assert(ipr_base==npr[is]);
    } // if s->non_local()
  }
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::update_twnl(void)
{
  // update arrays twnl[is][ipr][ig], dtwnl[is][ipr][j][ig],
  // following a change of cell dimensions
  // It is assumed that basis_ has been updated
  // It is assumed that nsp, npr[is], nquad[is] did not change since init

  tmap["update_twnl"].start();

  if ( nspnl == 0 ) return;

  const int ngwl = basis_.localsize();
  const double pi = M_PI;
  const double fpi = 4.0 * pi;
  const double s14pi = sqrt(1.0/fpi);
  const double s34pi = sqrt(3.0/fpi);
  const double s54pi = sqrt(5.0/fpi);
  const double s74pi = sqrt(7.0/fpi);
  const double s20pi = sqrt(20.0*pi);
  const double s20pi3 = sqrt(20.0*pi/3.0);
  const double s3 = sqrt(3.0);
  const double s32 = sqrt(1.5);
  const double s52 = sqrt(2.5);
  const double s15 = sqrt(15.0);

  const double *kpg   = basis_.kpg_ptr();
  const double *kpgi  = basis_.kpgi_ptr();
  const double *kpg_x = basis_.kpgx_ptr(0);
  const double *kpg_y = basis_.kpgx_ptr(1);
  const double *kpg_z = basis_.kpgx_ptr(2);

  // compute twnl and dtwnl
  for ( int is = 0; is < nsp; is++ )
  {
    Species *s = atoms_.species_list[is];

    int ilm = 0;
    for ( int iop = 0; iop < nop[is]; iop++ )
    {
      int l = s->l(iop);
      if ( l != lloc[is] )
      {
        if ( l == 0 )
        {
          if ( nquad[is] == 0 )
          {
            // Kleinman-Bylander

            // twnl[is][ipr][ig]
            const int ipr = ilm;
            // index = ig + ngwl*ipr, i.e. index = ig
            double *t0 = &twnl[is][ngwl*ipr];

            // dtwnl[is][ipr][ij][ngwl]
            // index = ig + ngwl * ( ij + 6 * ipr ), ipr = 0
            // i.e. index = ig + ij * ngwl
            double *dt0_xx = &dtwnl[is][ngwl*(0+6*ipr)];
            double *dt0_yy = &dtwnl[is][ngwl*(1+6*ipr)];
            double *dt0_zz = &dtwnl[is][ngwl*(2+6*ipr)];
            double *dt0_xy = &dtwnl[is][ngwl*(3+6*ipr)];
            double *dt0_yz = &dtwnl[is][ngwl*(4+6*ipr)];
            double *dt0_xz = &dtwnl[is][ngwl*(5+6*ipr)];
            // Special case k=G=0 is ok since kpgi[0] = 0.0 at k=G=0
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              s->dvnlg(iop,kpg[ig],v,dv);

              t0[ig] = s14pi * v;

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];

              const double tmp = kpgi[ig] * s14pi * dv;
              dt0_xx[ig] = tmp * tgx * tgx;
              dt0_yy[ig] = tmp * tgy * tgy;
              dt0_zz[ig] = tmp * tgz * tgz;
              dt0_xy[ig] = tmp * tgx * tgy;
              dt0_yz[ig] = tmp * tgy * tgz;
              dt0_xz[ig] = tmp * tgx * tgz;
            }
          }
          else
          {
            // semi-local
            for ( int iquad = 0; iquad < nquad[is]; iquad++ )
            {
              // twnl[is][ipr][ig]
              // ipr = iquad + nquad[is]*ilm, where ilm=0
              //     = iquad
              // index = ig + ngwl*iquad
              double *t0 = &twnl[is][ngwl*iquad];
              // dtwnl[is][ipr][j][ngwl]
              // index = ig + ngwl * ( ij + 6 * iquad)
              double *dt0_xx = &dtwnl[is][ngwl*(0+6*iquad)];
              double *dt0_yy = &dtwnl[is][ngwl*(1+6*iquad)];
              double *dt0_zz = &dtwnl[is][ngwl*(2+6*iquad)];
              double *dt0_xy = &dtwnl[is][ngwl*(3+6*iquad)];
              double *dt0_yz = &dtwnl[is][ngwl*(4+6*iquad)];
              double *dt0_xz = &dtwnl[is][ngwl*(5+6*iquad)];
              const double r = rquad[is][iquad];

              // G=0 element must be treated separately if k=0
              if ( basis_.real() && ctxt_.myrow() == 0)
              {
                // compute G=0 element separately to get
                // sin(Gr)/(Gr) limit -> 1
                const int ig = 0;
                // this task holds the G=0 element
                // I(l=0) = 4 pi j_l(G r) r
                // twnl[is][ipr][l][ig] = 4 pi j_0(Gr_i) r_i Ylm
                // j_0(Gr) = sin(Gr) / (Gr) == 1.0
                // Ylm = s14pi

                // next line: sin(Gr)/(Gr) replaced by 1
                t0[ig] = fpi * s14pi * r;

                // dtwnl = fpi s14pi G_i G_j / G (r cos(Gr)/G -sin(Gr)/G^2)
                // dtwnl = 0.0;
                dt0_xx[ig] = 0.0;
                dt0_yy[ig] = 0.0;
                dt0_zz[ig] = 0.0;
                dt0_xy[ig] = 0.0;
                dt0_yz[ig] = 0.0;
                dt0_xz[ig] = 0.0;
              }
              else
              {
                // Use the normal procedure
                const int ig = 0;
                // I(l=0) = 4 pi j_l(G r) r
                // twnl[is][ipr][l][ig] = 4 pi j_0(Gr_i) r_i Ylm
                // j_0(Gr) * r = sin(Gr) / G
                // Ylm = s14pi
                const double arg = kpg[ig] * r;
                // Note: for G=0, gi[0] = 0

                const double tgx = kpg_x[ig];
                const double tgy = kpg_y[ig];
                const double tgz = kpg_z[ig];
                const double tgi = kpgi[ig];
                const double tgi2 = tgi * tgi;

                const double ts = sin(arg);
                const double tc = cos(arg);

                t0[ig] = fpi * s14pi * ts * tgi;

                // dtwnl = fpi s14pi k+G_i k+G_j / k+G
                // (r cos(k+Gr)/k+G -sin(k+Gr)/k+G^2)
                const double tmp = fpi * s14pi * tgi2 * (r*tc - ts*tgi);
                dt0_xx[ig] = tmp * tgx * tgx;
                dt0_yy[ig] = tmp * tgy * tgy;
                dt0_zz[ig] = tmp * tgz * tgz;
                dt0_xy[ig] = tmp * tgx * tgy;
                dt0_yz[ig] = tmp * tgy * tgz;
                dt0_xz[ig] = tmp * tgx * tgz;
              }

              for ( int ig = 1; ig < ngwl; ig++ )
              {
                // I(l=0) = 4 pi j_l(|k+G| r) r
                // twnl[is][ipr][l][ig] = 4 pi j_0(|k+G|r_i) r_i Ylm
                // j_0(|k+G|r) * r = sin(|k+G|r) / |k+G|
                // Ylm = s14pi
                const double arg = kpg[ig] * r;

                const double tgx = kpg_x[ig];
                const double tgy = kpg_y[ig];
                const double tgz = kpg_z[ig];
                const double tgi = kpgi[ig];
                const double tgi2 = tgi * tgi;

                const double ts = sin(arg);
                const double tc = cos(arg);

                t0[ig] = fpi * s14pi * ts * tgi;

                // dtwnl = fpi s14pi (k+G)_i (k+G)_j /
                // |k+G| (r cos(|k+G|r)/|k+G| -sin(|k+G|r)/|k+G|^2)
                const double tmp = fpi * s14pi * tgi2 * (r*tc - ts*tgi);
                dt0_xx[ig] = tmp * tgx * tgx;
                dt0_yy[ig] = tmp * tgy * tgy;
                dt0_zz[ig] = tmp * tgz * tgz;
                dt0_xy[ig] = tmp * tgx * tgy;
                dt0_yz[ig] = tmp * tgy * tgz;
                dt0_xz[ig] = tmp * tgx * tgz;
              }
            } // for iquad
          }
          ilm += 2*l+1;
        }
        else if ( l == 1 )
        {
          if ( nquad[is] == 0 )
          {
            // Kleinman-Bylander

            // twnl[is][ipr][ig]
            // ipr = ilm
            const int ipr1 = ilm;
            const int ipr2 = ilm+1;
            const int ipr3 = ilm+2;
            // index = ig + ngwl*ilm
            double *t1 = &twnl[is][ngwl*ipr1];
            double *t2 = &twnl[is][ngwl*ipr2];
            double *t3 = &twnl[is][ngwl*ipr3];

            // dtwnl[is][ipr][ij][ngwl]
            // index = ig + ngwl * ( ij + 6 * ipr )
            double *dt1_xx = &dtwnl[is][ngwl*(0+6*ipr1)];
            double *dt1_yy = &dtwnl[is][ngwl*(1+6*ipr1)];
            double *dt1_zz = &dtwnl[is][ngwl*(2+6*ipr1)];
            double *dt1_xy = &dtwnl[is][ngwl*(3+6*ipr1)];
            double *dt1_yz = &dtwnl[is][ngwl*(4+6*ipr1)];
            double *dt1_xz = &dtwnl[is][ngwl*(5+6*ipr1)];

            double *dt2_xx = &dtwnl[is][ngwl*(0+6*ipr2)];
            double *dt2_yy = &dtwnl[is][ngwl*(1+6*ipr2)];
            double *dt2_zz = &dtwnl[is][ngwl*(2+6*ipr2)];
            double *dt2_xy = &dtwnl[is][ngwl*(3+6*ipr2)];
            double *dt2_yz = &dtwnl[is][ngwl*(4+6*ipr2)];
            double *dt2_xz = &dtwnl[is][ngwl*(5+6*ipr2)];

            double *dt3_xx = &dtwnl[is][ngwl*(0+6*ipr3)];
            double *dt3_yy = &dtwnl[is][ngwl*(1+6*ipr3)];
            double *dt3_zz = &dtwnl[is][ngwl*(2+6*ipr3)];
            double *dt3_xy = &dtwnl[is][ngwl*(3+6*ipr3)];
            double *dt3_yz = &dtwnl[is][ngwl*(4+6*ipr3)];
            double *dt3_xz = &dtwnl[is][ngwl*(5+6*ipr3)];

            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              const double tg = kpg[ig];
              s->dvnlg(iop,tg,v,dv);

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;

              const double tgi = kpgi[ig];
              const double tgi2 = tgi * tgi;

              const double y1 = s34pi * tgx * tgi;
              const double y2 = s34pi * tgy * tgi;
              const double y3 = s34pi * tgz * tgi;

              t1[ig]  = y1 * v;
              t2[ig]  = y2 * v;
              t3[ig]  = y3 * v;

              const double fac1 = - y1 * ( v - tg * dv ) * tgi2;
              // m=x
              dt1_xx[ig] = fac1 * tgx2 + v * y1;
              dt1_yy[ig] = fac1 * tgy2;
              dt1_zz[ig] = fac1 * tgz2;
              dt1_xy[ig] = fac1 * tgx * tgy;
              dt1_yz[ig] = fac1 * tgy * tgz;
              dt1_xz[ig] = fac1 * tgx * tgz;

              const double fac2 = - y2 * ( v - tg * dv ) * tgi2;
              // m=y
              dt2_xx[ig] = fac2 * tgx2;
              dt2_yy[ig] = fac2 * tgy2 + v * y2;
              dt2_zz[ig] = fac2 * tgz2;
              dt2_xy[ig] = fac2 * tgx * tgy + v * y1;
              dt2_yz[ig] = fac2 * tgy * tgz;
              dt2_xz[ig] = fac2 * tgx * tgz;

              const double fac3 = - y3 * ( v - tg * dv ) * tgi2;
              // m=z
              dt3_xx[ig] = fac3 * tgx2;
              dt3_yy[ig] = fac3 * tgy2;
              dt3_zz[ig] = fac3 * tgz2 + v * y3;
              dt3_xy[ig] = fac3 * tgx * tgy;
              dt3_yz[ig] = fac3 * tgy * tgz + v * y2;
              dt3_xz[ig] = fac3 * tgx * tgz + v * y1;
            }
          }
          else
          {
            // semi-local
            for ( int iquad = 0; iquad < nquad[is]; iquad++ )
            {
              // twnl[is][ipr][ig]
              // index = ig + ngwl*(iquad+nquad[is]*ilm)
              const int ipr1 = iquad+nquad[is]*ilm;
              const int ipr2 = iquad+nquad[is]*(ilm+1);
              const int ipr3 = iquad+nquad[is]*(ilm+2);
              double *t1 = &twnl[is][ngwl*ipr1];
              double *t2 = &twnl[is][ngwl*ipr2];
              double *t3 = &twnl[is][ngwl*ipr3];

              // dtwnl[is][ipr][j][ngwl]
              // index = ig + ngwl * ( ij + 6 * ipr )
              double *dt1_xx = &dtwnl[is][ngwl*(0+6*ipr1)];
              double *dt1_yy = &dtwnl[is][ngwl*(1+6*ipr1)];
              double *dt1_zz = &dtwnl[is][ngwl*(2+6*ipr1)];
              double *dt1_xy = &dtwnl[is][ngwl*(3+6*ipr1)];
              double *dt1_yz = &dtwnl[is][ngwl*(4+6*ipr1)];
              double *dt1_xz = &dtwnl[is][ngwl*(5+6*ipr1)];

              double *dt2_xx = &dtwnl[is][ngwl*(0+6*ipr2)];
              double *dt2_yy = &dtwnl[is][ngwl*(1+6*ipr2)];
              double *dt2_zz = &dtwnl[is][ngwl*(2+6*ipr2)];
              double *dt2_xy = &dtwnl[is][ngwl*(3+6*ipr2)];
              double *dt2_yz = &dtwnl[is][ngwl*(4+6*ipr2)];
              double *dt2_xz = &dtwnl[is][ngwl*(5+6*ipr2)];

              double *dt3_xx = &dtwnl[is][ngwl*(0+6*ipr3)];
              double *dt3_yy = &dtwnl[is][ngwl*(1+6*ipr3)];
              double *dt3_zz = &dtwnl[is][ngwl*(2+6*ipr3)];
              double *dt3_xy = &dtwnl[is][ngwl*(3+6*ipr3)];
              double *dt3_yz = &dtwnl[is][ngwl*(4+6*ipr3)];
              double *dt3_xz = &dtwnl[is][ngwl*(5+6*ipr3)];

              const double r = rquad[is][iquad];
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                double v = 0.0, dv = 0.0;
                // j_1(Gr) = (1/(Gr))*(sin(Gr)/(Gr)-cos(Gr))
                const double tg = kpg[ig];
                const double z = tg * r;
                if ( z != 0.0 )
                {
                  const double zi = 1.0 / z;
                  const double c = cos(z);
                  const double s = sin(z);
                  const double j1 = ( s * zi - c ) * zi;
                  const double dj1 =
                    ( 2.0 * z * c + ( z*z - 2.0 ) * s ) * zi*zi*zi;
                  // v = 4 pi j1(Gr) r
                  v = fpi * j1 * r;
                  // dv = d/dG v = 4 pi dj1(Gr)/d(Gr) d(Gr)/dG r
                  //    = 4 pi dj1 r^2
                  dv = fpi * dj1 * r * r;
                }

                const double tgx = kpg_x[ig];
                const double tgy = kpg_y[ig];
                const double tgz = kpg_z[ig];
                const double tgx2 = tgx * tgx;
                const double tgy2 = tgy * tgy;
                const double tgz2 = tgz * tgz;

                const double tgi = kpgi[ig];
                const double tgi2 = tgi * tgi;

                const double y1 = s34pi * tgx * tgi;
                const double y2 = s34pi * tgy * tgi;
                const double y3 = s34pi * tgz * tgi;

                t1[ig]  = y1 * v;
                t2[ig]  = y2 * v;
                t3[ig]  = y3 * v;

                const double fac1 = - y1 * ( v - tg * dv ) * tgi2;
                // m=x
                dt1_xx[ig] = fac1 * tgx2 + v * y1;
                dt1_yy[ig] = fac1 * tgy2;
                dt1_zz[ig] = fac1 * tgz2;
                dt1_xy[ig] = fac1 * tgx * tgy;
                dt1_yz[ig] = fac1 * tgy * tgz;
                dt1_xz[ig] = fac1 * tgx * tgz;

                const double fac2 = - y2 * ( v - tg * dv ) * tgi2;
                // m=y
                dt2_xx[ig] = fac2 * tgx2;
                dt2_yy[ig] = fac2 * tgy2 + v * y2;
                dt2_zz[ig] = fac2 * tgz2;
                dt2_xy[ig] = fac2 * tgx * tgy + v * y1;
                dt2_yz[ig] = fac2 * tgy * tgz;
                dt2_xz[ig] = fac2 * tgx * tgz;

                const double fac3 = - y3 * ( v - tg * dv ) * tgi2;
                // m=z
                dt3_xx[ig] = fac3 * tgx2;
                dt3_yy[ig] = fac3 * tgy2;
                dt3_zz[ig] = fac3 * tgz2 + v * y3;
                dt3_xy[ig] = fac3 * tgx * tgy;
                dt3_yz[ig] = fac3 * tgy * tgz + v * y2;
                dt3_xz[ig] = fac3 * tgx * tgz + v * y1;
              } // ig
            }
          }
          ilm += 2*l+1;
        }
        else if ( l == 2 )
        {
          if ( nquad[is] == 0 )
          {
            // Kleinman-Bylander
            const int ipr4 = ilm;
            const int ipr5 = ilm+1;
            const int ipr6 = ilm+2;
            const int ipr7 = ilm+3;
            const int ipr8 = ilm+4;

            double *t4 = &twnl[is][ngwl*ipr4];
            double *t5 = &twnl[is][ngwl*ipr5];
            double *t6 = &twnl[is][ngwl*ipr6];
            double *t7 = &twnl[is][ngwl*ipr7];
            double *t8 = &twnl[is][ngwl*ipr8];

            // dtwnl[is][ipr][ij][ngwl]
            // index = ig + ngwl * ( ij + 6 * ipr )
            double *dt4_xx = &dtwnl[is][ngwl*(0+6*ipr4)];
            double *dt4_yy = &dtwnl[is][ngwl*(1+6*ipr4)];
            double *dt4_zz = &dtwnl[is][ngwl*(2+6*ipr4)];
            double *dt4_xy = &dtwnl[is][ngwl*(3+6*ipr4)];
            double *dt4_yz = &dtwnl[is][ngwl*(4+6*ipr4)];
            double *dt4_xz = &dtwnl[is][ngwl*(5+6*ipr4)];

            double *dt5_xx = &dtwnl[is][ngwl*(0+6*ipr5)];
            double *dt5_yy = &dtwnl[is][ngwl*(1+6*ipr5)];
            double *dt5_zz = &dtwnl[is][ngwl*(2+6*ipr5)];
            double *dt5_xy = &dtwnl[is][ngwl*(3+6*ipr5)];
            double *dt5_yz = &dtwnl[is][ngwl*(4+6*ipr5)];
            double *dt5_xz = &dtwnl[is][ngwl*(5+6*ipr5)];

            double *dt6_xx = &dtwnl[is][ngwl*(0+6*ipr6)];
            double *dt6_yy = &dtwnl[is][ngwl*(1+6*ipr6)];
            double *dt6_zz = &dtwnl[is][ngwl*(2+6*ipr6)];
            double *dt6_xy = &dtwnl[is][ngwl*(3+6*ipr6)];
            double *dt6_yz = &dtwnl[is][ngwl*(4+6*ipr6)];
            double *dt6_xz = &dtwnl[is][ngwl*(5+6*ipr6)];

            double *dt7_xx = &dtwnl[is][ngwl*(0+6*ipr7)];
            double *dt7_yy = &dtwnl[is][ngwl*(1+6*ipr7)];
            double *dt7_zz = &dtwnl[is][ngwl*(2+6*ipr7)];
            double *dt7_xy = &dtwnl[is][ngwl*(3+6*ipr7)];
            double *dt7_yz = &dtwnl[is][ngwl*(4+6*ipr7)];
            double *dt7_xz = &dtwnl[is][ngwl*(5+6*ipr7)];

            double *dt8_xx = &dtwnl[is][ngwl*(0+6*ipr8)];
            double *dt8_yy = &dtwnl[is][ngwl*(1+6*ipr8)];
            double *dt8_zz = &dtwnl[is][ngwl*(2+6*ipr8)];
            double *dt8_xy = &dtwnl[is][ngwl*(3+6*ipr8)];
            double *dt8_yz = &dtwnl[is][ngwl*(4+6*ipr8)];
            double *dt8_xz = &dtwnl[is][ngwl*(5+6*ipr8)];

            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              const double tg = kpg[ig];

              s->dvnlg(iop,tg,v,dv);

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;

              const double tgi = kpgi[ig];
              const double tgi2 = tgi * tgi;

              const double tgxx = tgx2 * tgi2;
              const double tgyy = tgy2 * tgi2;
              const double tgzz = tgz2 * tgi2;
              const double tgxy = tgx * tgy * tgi2;
              const double tgyz = tgy * tgz * tgi2;
              const double tgxz = tgx * tgz * tgi2;

              const double y4 = s54pi * 0.5 * (3.0 * tgzz - 1.0 );
              const double y5 = s54pi * 0.5 * s3 * ( tgxx - tgyy );
              const double y6 = s54pi * s3 * tgxy;
              const double y7 = s54pi * s3 * tgyz;
              const double y8 = s54pi * s3 * tgxz;

              const double y1x = s34pi * tgx * tgi;
              const double y1y = s34pi * tgy * tgi;
              const double y1z = s34pi * tgz * tgi;

              const double dx_xx = y1x * tgxx - y1x;
              const double dx_yy = y1x * tgyy;
              const double dx_zz = y1x * tgzz;
              const double dx_xy = y1x * tgxy;
              const double dx_yz = y1x * tgyz;
              const double dx_xz = y1x * tgxz;

              const double dy_xx = y1y * tgxx;
              const double dy_yy = y1y * tgyy - y1y;
              const double dy_zz = y1y * tgzz;
              const double dy_xy = y1y * tgxy - y1x;
              const double dy_yz = y1y * tgyz;
              const double dy_xz = y1y * tgxz;

              const double dz_xx = y1z * tgxx;
              const double dz_yy = y1z * tgyy;
              const double dz_zz = y1z * tgzz - y1z;
              const double dz_xy = y1z * tgxy;
              const double dz_yz = y1z * tgyz - y1y;
              const double dz_xz = y1z * tgxz - y1x;

              t4[ig]  = y4 * v;
              t5[ig]  = y5 * v;
              t6[ig]  = y6 * v;
              t7[ig]  = y7 * v;
              t8[ig]  = y8 * v;

              // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
              dt4_xx[ig] = -(v * s20pi * dz_xx * y1z - y4 * dv * tg * tgxx);
              dt4_yy[ig] = -(v * s20pi * dz_yy * y1z - y4 * dv * tg * tgyy);
              dt4_zz[ig] = -(v * s20pi * dz_zz * y1z - y4 * dv * tg * tgzz);
              dt4_xy[ig] = -(v * s20pi * dz_xy * y1z - y4 * dv * tg * tgxy);
              dt4_yz[ig] = -(v * s20pi * dz_yz * y1z - y4 * dv * tg * tgyz);
              dt4_xz[ig] = -(v * s20pi * dz_xz * y1z - y4 * dv * tg * tgxz);

              // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
              dt5_xx[ig] = -(v * s20pi3 * (y1x * dx_xx - y1y * dy_xx) - y5 * dv * tg * tgxx);
              dt5_yy[ig] = -(v * s20pi3 * (y1x * dx_yy - y1y * dy_yy) - y5 * dv * tg * tgyy);
              dt5_zz[ig] = -(v * s20pi3 * (y1x * dx_zz - y1y * dy_zz) - y5 * dv * tg * tgzz);
              dt5_xy[ig] = -(v * s20pi3 * (y1x * dx_xy - y1y * dy_xy) - y5 * dv * tg * tgxy);
              dt5_yz[ig] = -(v * s20pi3 * (y1x * dx_yz - y1y * dy_yz) - y5 * dv * tg * tgyz);
              dt5_xz[ig] = -(v * s20pi3 * (y1x * dx_xz - y1y * dy_xz) - y5 * dv * tg * tgxz);

              // y6 = s54pi sqrt(3) x y / r^2
              dt6_xx[ig] = -(v * s20pi3 * (dx_xx * y1y + y1x * dy_xx) - y6 * dv * tg * tgxx);
              dt6_yy[ig] = -(v * s20pi3 * (dx_yy * y1y + y1x * dy_yy) - y6 * dv * tg * tgyy);
              dt6_zz[ig] = -(v * s20pi3 * (dx_zz * y1y + y1x * dy_zz) - y6 * dv * tg * tgzz);
              dt6_xy[ig] = -(v * s20pi3 * (dx_xy * y1y + y1x * dy_xy) - y6 * dv * tg * tgxy);
              dt6_yz[ig] = -(v * s20pi3 * (dx_yz * y1y + y1x * dy_yz) - y6 * dv * tg * tgyz);
              dt6_xz[ig] = -(v * s20pi3 * (dx_xz * y1y + y1x * dy_xz) - y6 * dv * tg * tgxz);

              // y7 = s54pi sqrt(3) y z / r^2
              dt7_xx[ig] = -(v * s20pi3 * (dy_xx * y1z + y1y * dz_xx) - y7 * dv * tg * tgxx);
              dt7_yy[ig] = -(v * s20pi3 * (dy_yy * y1z + y1y * dz_yy) - y7 * dv * tg * tgyy);
              dt7_zz[ig] = -(v * s20pi3 * (dy_zz * y1z + y1y * dz_zz) - y7 * dv * tg * tgzz);
              dt7_xy[ig] = -(v * s20pi3 * (dy_xy * y1z + y1y * dz_xy) - y7 * dv * tg * tgxy);
              dt7_yz[ig] = -(v * s20pi3 * (dy_yz * y1z + y1y * dz_yz) - y7 * dv * tg * tgyz);
              dt7_xz[ig] = -(v * s20pi3 * (dy_xz * y1z + y1y * dz_xz) - y7 * dv * tg * tgxz);

              // y8 = s54pi sqrt(3) z x / r^2
              dt8_xx[ig] = -(v * s20pi3 * (dx_xx * y1z + y1x * dz_xx) - y8 * dv * tg * tgxx);
              dt8_yy[ig] = -(v * s20pi3 * (dx_yy * y1z + y1x * dz_yy) - y8 * dv * tg * tgyy);
              dt8_zz[ig] = -(v * s20pi3 * (dx_zz * y1z + y1x * dz_zz) - y8 * dv * tg * tgzz);
              dt8_xy[ig] = -(v * s20pi3 * (dx_xy * y1z + y1x * dz_xy) - y8 * dv * tg * tgxy);
              dt8_yz[ig] = -(v * s20pi3 * (dx_yz * y1z + y1x * dz_yz) - y8 * dv * tg * tgyz);
              dt8_xz[ig] = -(v * s20pi3 * (dx_xz * y1z + y1x * dz_xz) - y8 * dv * tg * tgxz);
            }
          }
          else
          {
            // semi-local
            for ( int iquad = 0; iquad < nquad[is]; iquad++ )
            {
              // twnl[is][ipr][ig]
              // ipr = iquad+nquad[is]*ilm
              // index = ig + ngwl*ipr
              const int ipr4 = iquad+nquad[is]*ilm;
              const int ipr5 = iquad+nquad[is]*(ilm+1);
              const int ipr6 = iquad+nquad[is]*(ilm+2);
              const int ipr7 = iquad+nquad[is]*(ilm+3);
              const int ipr8 = iquad+nquad[is]*(ilm+4);

              double *t4 = &twnl[is][ngwl*ipr4];
              double *t5 = &twnl[is][ngwl*ipr5];
              double *t6 = &twnl[is][ngwl*ipr6];
              double *t7 = &twnl[is][ngwl*ipr7];
              double *t8 = &twnl[is][ngwl*ipr8];

              // dtwnl[is][ipr][ij][ngwl]
              // index = ig + ngwl * ( ij + 6 * ipr )
              double *dt4_xx = &dtwnl[is][ngwl*(0+6*ipr4)];
              double *dt4_yy = &dtwnl[is][ngwl*(1+6*ipr4)];
              double *dt4_zz = &dtwnl[is][ngwl*(2+6*ipr4)];
              double *dt4_xy = &dtwnl[is][ngwl*(3+6*ipr4)];
              double *dt4_yz = &dtwnl[is][ngwl*(4+6*ipr4)];
              double *dt4_xz = &dtwnl[is][ngwl*(5+6*ipr4)];

              double *dt5_xx = &dtwnl[is][ngwl*(0+6*ipr5)];
              double *dt5_yy = &dtwnl[is][ngwl*(1+6*ipr5)];
              double *dt5_zz = &dtwnl[is][ngwl*(2+6*ipr5)];
              double *dt5_xy = &dtwnl[is][ngwl*(3+6*ipr5)];
              double *dt5_yz = &dtwnl[is][ngwl*(4+6*ipr5)];
              double *dt5_xz = &dtwnl[is][ngwl*(5+6*ipr5)];

              double *dt6_xx = &dtwnl[is][ngwl*(0+6*ipr6)];
              double *dt6_yy = &dtwnl[is][ngwl*(1+6*ipr6)];
              double *dt6_zz = &dtwnl[is][ngwl*(2+6*ipr6)];
              double *dt6_xy = &dtwnl[is][ngwl*(3+6*ipr6)];
              double *dt6_yz = &dtwnl[is][ngwl*(4+6*ipr6)];
              double *dt6_xz = &dtwnl[is][ngwl*(5+6*ipr6)];

              double *dt7_xx = &dtwnl[is][ngwl*(0+6*ipr7)];
              double *dt7_yy = &dtwnl[is][ngwl*(1+6*ipr7)];
              double *dt7_zz = &dtwnl[is][ngwl*(2+6*ipr7)];
              double *dt7_xy = &dtwnl[is][ngwl*(3+6*ipr7)];
              double *dt7_yz = &dtwnl[is][ngwl*(4+6*ipr7)];
              double *dt7_xz = &dtwnl[is][ngwl*(5+6*ipr7)];

              double *dt8_xx = &dtwnl[is][ngwl*(0+6*ipr8)];
              double *dt8_yy = &dtwnl[is][ngwl*(1+6*ipr8)];
              double *dt8_zz = &dtwnl[is][ngwl*(2+6*ipr8)];
              double *dt8_xy = &dtwnl[is][ngwl*(3+6*ipr8)];
              double *dt8_yz = &dtwnl[is][ngwl*(4+6*ipr8)];
              double *dt8_xz = &dtwnl[is][ngwl*(5+6*ipr8)];

              const double r = rquad[is][iquad];
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                double v = 0.0, dv = 0.0;
                // j_2(z) = (3/z^3-1/z) sin(z) - 3/z^2 cos(z)
                // j_2(z) = (1/z)*((3/z^2-1)*sin(z) - (3/z) cos(z) )

                const double tg = kpg[ig];
                const double z = tg * r;
                if ( z != 0.0 )
                {
                  const double zi = 1.0 / z;
                  const double c = cos(z);
                  const double s = sin(z);
                  const double j2 = ((3.0*zi*zi - 1.0) * s - 3.0*zi * c ) * zi;
                  const double z2 = z * z;
                  const double dj2 =
                    ( (4.0 * z2 - 9.0) * s + z*(9.0-z2) * c ) / (z2*z2) ;
                  // v = 4 pi j2(Gr) r
                  v = fpi * j2 * r;
                  // dv = d/dG v = 4 pi dj2(Gr)/d(Gr) d(Gr)/dG r
                  //    = 4 pi dj2 r^2
                  dv = fpi * dj2 * r * r;
                }

                const double tgx = kpg_x[ig];
                const double tgy = kpg_y[ig];
                const double tgz = kpg_z[ig];
                const double tgx2 = tgx * tgx;
                const double tgy2 = tgy * tgy;
                const double tgz2 = tgz * tgz;

                const double tgi = kpgi[ig];
                const double tgi2 = tgi * tgi;

                const double tgxx = tgx2 * tgi2;
                const double tgyy = tgy2 * tgi2;
                const double tgzz = tgz2 * tgi2;
                const double tgxy = tgx * tgy * tgi2;
                const double tgyz = tgy * tgz * tgi2;
                const double tgxz = tgx * tgz * tgi2;

                const double y4 = s54pi * 0.5 * (3.0 * tgzz - 1.0 );
                const double y5 = s54pi * 0.5 * s3 * ( tgxx - tgyy );
                const double y6 = s54pi * s3 * tgxy;
                const double y7 = s54pi * s3 * tgyz;
                const double y8 = s54pi * s3 * tgxz;

                const double y1x = s34pi * tgx * tgi;
                const double y1y = s34pi * tgy * tgi;
                const double y1z = s34pi * tgz * tgi;

                const double dx_xx = y1x * tgxx - y1x;
                const double dx_yy = y1x * tgyy;
                const double dx_zz = y1x * tgzz;
                const double dx_xy = y1x * tgxy;
                const double dx_yz = y1x * tgyz;
                const double dx_xz = y1x * tgxz;

                const double dy_xx = y1y * tgxx;
                const double dy_yy = y1y * tgyy - y1y;
                const double dy_zz = y1y * tgzz;
                const double dy_xy = y1y * tgxy - y1x;
                const double dy_yz = y1y * tgyz;
                const double dy_xz = y1y * tgxz;

                const double dz_xx = y1z * tgxx;
                const double dz_yy = y1z * tgyy;
                const double dz_zz = y1z * tgzz - y1z;
                const double dz_xy = y1z * tgxy;
                const double dz_yz = y1z * tgyz - y1y;
                const double dz_xz = y1z * tgxz - y1x;

                t4[ig]  = y4 * v;
                t5[ig]  = y5 * v;
                t6[ig]  = y6 * v;
                t7[ig]  = y7 * v;
                t8[ig]  = y8 * v;

                // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
                dt4_xx[ig] = -(v * s20pi * dz_xx * y1z - y4 * dv * tg * tgxx);
                dt4_yy[ig] = -(v * s20pi * dz_yy * y1z - y4 * dv * tg * tgyy);
                dt4_zz[ig] = -(v * s20pi * dz_zz * y1z - y4 * dv * tg * tgzz);
                dt4_xy[ig] = -(v * s20pi * dz_xy * y1z - y4 * dv * tg * tgxy);
                dt4_yz[ig] = -(v * s20pi * dz_yz * y1z - y4 * dv * tg * tgyz);
                dt4_xz[ig] = -(v * s20pi * dz_xz * y1z - y4 * dv * tg * tgxz);

                // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
                dt5_xx[ig] = -(v * s20pi3 * (y1x * dx_xx - y1y * dy_xx) - y5 * dv * tg * tgxx);
                dt5_yy[ig] = -(v * s20pi3 * (y1x * dx_yy - y1y * dy_yy) - y5 * dv * tg * tgyy);
                dt5_zz[ig] = -(v * s20pi3 * (y1x * dx_zz - y1y * dy_zz) - y5 * dv * tg * tgzz);
                dt5_xy[ig] = -(v * s20pi3 * (y1x * dx_xy - y1y * dy_xy) - y5 * dv * tg * tgxy);
                dt5_yz[ig] = -(v * s20pi3 * (y1x * dx_yz - y1y * dy_yz) - y5 * dv * tg * tgyz);
                dt5_xz[ig] = -(v * s20pi3 * (y1x * dx_xz - y1y * dy_xz) - y5 * dv * tg * tgxz);

                // y6 = s54pi sqrt(3) x y / r^2
                dt6_xx[ig] = -(v * s20pi3 * (dx_xx * y1y + y1x * dy_xx) - y6 * dv * tg * tgxx);
                dt6_yy[ig] = -(v * s20pi3 * (dx_yy * y1y + y1x * dy_yy) - y6 * dv * tg * tgyy);
                dt6_zz[ig] = -(v * s20pi3 * (dx_zz * y1y + y1x * dy_zz) - y6 * dv * tg * tgzz);
                dt6_xy[ig] = -(v * s20pi3 * (dx_xy * y1y + y1x * dy_xy) - y6 * dv * tg * tgxy);
                dt6_yz[ig] = -(v * s20pi3 * (dx_yz * y1y + y1x * dy_yz) - y6 * dv * tg * tgyz);
                dt6_xz[ig] = -(v * s20pi3 * (dx_xz * y1y + y1x * dy_xz) - y6 * dv * tg * tgxz);

                // y7 = s54pi sqrt(3) y z / r^2
                dt7_xx[ig] = -(v * s20pi3 * (dy_xx * y1z + y1y * dz_xx) - y7 * dv * tg * tgxx);
                dt7_yy[ig] = -(v * s20pi3 * (dy_yy * y1z + y1y * dz_yy) - y7 * dv * tg * tgyy);
                dt7_zz[ig] = -(v * s20pi3 * (dy_zz * y1z + y1y * dz_zz) - y7 * dv * tg * tgzz);
                dt7_xy[ig] = -(v * s20pi3 * (dy_xy * y1z + y1y * dz_xy) - y7 * dv * tg * tgxy);
                dt7_yz[ig] = -(v * s20pi3 * (dy_yz * y1z + y1y * dz_yz) - y7 * dv * tg * tgyz);
                dt7_xz[ig] = -(v * s20pi3 * (dy_xz * y1z + y1y * dz_xz) - y7 * dv * tg * tgxz);

                // y8 = s54pi sqrt(3) z x / r^2
                dt8_xx[ig] = -(v * s20pi3 * (dx_xx * y1z + y1x * dz_xx) - y8 * dv * tg * tgxx);
                dt8_yy[ig] = -(v * s20pi3 * (dx_yy * y1z + y1x * dz_yy) - y8 * dv * tg * tgyy);
                dt8_zz[ig] = -(v * s20pi3 * (dx_zz * y1z + y1x * dz_zz) - y8 * dv * tg * tgzz);
                dt8_xy[ig] = -(v * s20pi3 * (dx_xy * y1z + y1x * dz_xy) - y8 * dv * tg * tgxy);
                dt8_yz[ig] = -(v * s20pi3 * (dx_yz * y1z + y1x * dz_yz) - y8 * dv * tg * tgyz);
                dt8_xz[ig] = -(v * s20pi3 * (dx_xz * y1z + y1x * dz_xz) - y8 * dv * tg * tgxz);
              } // ig
            } // iquad
          }
          ilm += 2*l+1;
        }
        else if ( l == 3 )
        {
          // only implemented for Kleiman-Bylander type
          assert(nquad[is] == 0);
          // Kleinman-Bylander
          const int ipr09 = ilm;
          const int ipr10 = ilm + 1;
          const int ipr11 = ilm + 2;
          const int ipr12 = ilm + 3;
          const int ipr13 = ilm + 4;
          const int ipr14 = ilm + 5;
          const int ipr15 = ilm + 6;

          double *t09 = &twnl[is][ngwl * ipr09];
          double *t10 = &twnl[is][ngwl * ipr10];
          double *t11 = &twnl[is][ngwl * ipr11];
          double *t12 = &twnl[is][ngwl * ipr12];
          double *t13 = &twnl[is][ngwl * ipr13];
          double *t14 = &twnl[is][ngwl * ipr14];
          double *t15 = &twnl[is][ngwl * ipr15];

          // dtwnl[is][ipr][ij][ngwl]
          // index = ig + ngwl * ( ij + 6 * ipr )
          double *dt09_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr09 )];
          double *dt09_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr09 )];
          double *dt09_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr09 )];
          double *dt09_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr09 )];
          double *dt09_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr09 )];
          double *dt09_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr09 )];

          double *dt10_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr10 )];
          double *dt10_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr10 )];
          double *dt10_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr10 )];
          double *dt10_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr10 )];
          double *dt10_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr10 )];
          double *dt10_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr10 )];

          double *dt11_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr11 )];
          double *dt11_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr11 )];
          double *dt11_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr11 )];
          double *dt11_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr11 )];
          double *dt11_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr11 )];
          double *dt11_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr11 )];

          double *dt12_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr12 )];
          double *dt12_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr12 )];
          double *dt12_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr12 )];
          double *dt12_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr12 )];
          double *dt12_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr12 )];
          double *dt12_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr12 )];

          double *dt13_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr13 )];
          double *dt13_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr13 )];
          double *dt13_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr13 )];
          double *dt13_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr13 )];
          double *dt13_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr13 )];
          double *dt13_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr13 )];

          double *dt14_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr14 )];
          double *dt14_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr14 )];
          double *dt14_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr14 )];
          double *dt14_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr14 )];
          double *dt14_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr14 )];
          double *dt14_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr14 )];

          double *dt15_xx = &dtwnl[is][ngwl * ( 0 + 6 * ipr15 )];
          double *dt15_yy = &dtwnl[is][ngwl * ( 1 + 6 * ipr15 )];
          double *dt15_zz = &dtwnl[is][ngwl * ( 2 + 6 * ipr15 )];
          double *dt15_xy = &dtwnl[is][ngwl * ( 3 + 6 * ipr15 )];
          double *dt15_yz = &dtwnl[is][ngwl * ( 4 + 6 * ipr15 )];
          double *dt15_xz = &dtwnl[is][ngwl * ( 5 + 6 * ipr15 )];

          for ( int ig = 0; ig < ngwl; ig++ )
          {
            double v, dv;
            const double tg = kpg[ig];

            s->dvnlg(iop,tg,v,dv);

            const double tgx = kpg_x[ig];
            const double tgy = kpg_y[ig];
            const double tgz = kpg_z[ig];
            const double tgx2 = tgx * tgx;
            const double tgy2 = tgy * tgy;
            const double tgz2 = tgz * tgz;
            const double tgx3 = tgx * tgx2;
            const double tgy3 = tgy * tgy2;
            const double tgz3 = tgz * tgz2;

            const double tgi = kpgi[ig];
            const double tgi2 = tgi * tgi;
            const double tgi3 = tgi * tgi2;

            const double tgxx = tgx2 * tgi2;
            const double tgyy = tgy2 * tgi2;
            const double tgzz = tgz2 * tgi2;
            const double tgxy = tgx * tgy * tgi2;
            const double tgyz = tgy * tgz * tgi2;
            const double tgxz = tgx * tgz * tgi2;

            const double tgxxx = tgx3 * tgi3;
            const double tgyyy = tgy3 * tgi3;
            const double tgzzz = tgz3 * tgi3;
            const double tgxyy = tgx * tgy2 * tgi3;
            const double tgxzz = tgx * tgz2 * tgi3;
            const double tgyxx = tgy * tgx2 * tgi3;
            const double tgyzz = tgy * tgz2 * tgi3;
            const double tgzxx = tgz * tgx2 * tgi3;
            const double tgzyy = tgz * tgy2 * tgi3;
            const double tgxyz = tgx * tgy * tgz * tgi3;

            const double factor = 0.5 * s74pi;

            const double y09 = factor * s52 * ( 3.0 * tgyxx - tgyyy );
            const double y10 = factor * s15 * 2.0 * tgxyz;
            const double y11 = factor * s32 * ( 4.0 * tgyzz - tgyxx - tgyyy );
            const double y12 = factor
              * ( 2.0 * tgzzz - 3.0 * ( tgzxx + tgzyy ) );
            const double y13 = factor * s32 * ( 4.0 * tgxzz - tgxxx - tgxyy );
            const double y14 = factor * s15 * ( tgzxx - tgzyy );
            const double y15 = factor * s52 * ( tgxxx - 3.0 * tgxyy );

            // derivative of x^3/r^3 w.r.t. x, y, and z
            const double dx_xxx = 3.0 * tgxx * ( 1.0 - tgxx ) * tgi;
            const double dy_xxx = -3.0 * tgxx * tgxy * tgi;
            const double dz_xxx = -3.0 * tgxx * tgxz * tgi;
            // derivative of y^3/r^3 w.r.t. x, y, and z
            const double dx_yyy = -3.0 * tgyy * tgxy * tgi;
            const double dy_yyy = 3.0 * tgyy * ( 1.0 - tgyy ) * tgi;
            const double dz_yyy = -3.0 * tgyy * tgyz * tgi;
            // derivative of z^3/r^3 w.r.t. x, y, and z
            const double dx_zzz = -3.0 * tgzz * tgxz * tgi;
            const double dy_zzz = -3.0 * tgzz * tgyz * tgi;
            const double dz_zzz = 3.0 * tgzz * ( 1.0 - tgzz ) * tgi;
            // derivative of xy^2/r^3 w.r.t. x, y, and z
            const double dx_xyy = tgyy * ( 1.0 - 3.0 * tgxx ) * tgi;
            const double dy_xyy = tgxy * ( 2.0 - 3.0 * tgyy ) * tgi;
            const double dz_xyy = -3.0 * tgyy * tgxz * tgi;
            // derivative of xz^2/r^3 w.r.t. x, y, and z
            const double dx_xzz = tgzz * ( 1.0 - 3.0 * tgxx ) * tgi;
            const double dy_xzz = -3.0 * tgzz * tgxy * tgi;
            const double dz_xzz = tgxz * ( 2.0 - 3.0 * tgzz ) * tgi;
            // derivative of yx^2/r^3 w.r.t. x, y, and z
            const double dx_yxx = tgxy * ( 2.0 - 3.0 * tgxx ) * tgi;
            const double dy_yxx = tgxx * ( 1.0 - 3.0 * tgyy ) * tgi;
            const double dz_yxx = -3.0 * tgxx * tgyz * tgi;
            // derivative of yz^2/r^3 w.r.t. x, y, and z
            const double dx_yzz = -3.0 * tgzz * tgxy * tgi;
            const double dy_yzz = tgzz * ( 1.0 - 3.0 * tgyy ) * tgi;
            const double dz_yzz = tgyz * ( 2.0 - 3.0 * tgzz ) * tgi;
            // derivative of zx^2/r^3 w.r.t. x, y, and z
            const double dx_zxx = tgxz * ( 2.0 - 3.0 * tgxx ) * tgi;
            const double dy_zxx = -3.0 * tgxx * tgyz * tgi;
            const double dz_zxx = tgxx * ( 1.0 - 3.0 * tgzz ) * tgi;
            // derivative of zy^2/r^3 w.r.t. x, y, and z
            const double dx_zyy = -3.0 * tgyy * tgxz * tgi;
            const double dy_zyy = tgyz * ( 2.0 - 3.0 * tgyy ) * tgi;
            const double dz_zyy = tgyy * ( 1.0 - 3.0 * tgzz ) * tgi;
            // derivative of xyz/r^3 w.r.t. x, y, and z
            const double dx_xyz = tgyz * ( 1 - 3.0 * tgxx ) * tgi;
            const double dy_xyz = tgxz * ( 1 - 3.0 * tgyy ) * tgi;
            const double dz_xyz = tgxy * ( 1 - 3.0 * tgzz ) * tgi;

            // derivatives of spherical harmonics

            // y9 = factor * s52 * ( 3.0 * tgyxx - tgyyy );
            const double dx_y09 = factor * s52 * ( 3.0 * dx_yxx - dx_yyy );
            const double dy_y09 = factor * s52 * ( 3.0 * dy_yxx - dy_yyy );
            const double dz_y09 = factor * s52 * ( 3.0 * dz_yxx - dz_yyy );
            // y10 = factor * s15 * 2.0 * tgxyz;
            const double dx_y10 = factor * s15 * 2.0 * dx_xyz;
            const double dy_y10 = factor * s15 * 2.0 * dy_xyz;
            const double dz_y10 = factor * s15 * 2.0 * dz_xyz;
            // y11 = factor * s32 * ( 4.0 * tgyzz - tgyxx - tgyyy );
            const double dx_y11 = factor * s32 * ( 4.0 * dx_yzz - dx_yxx - dx_yyy );
            const double dy_y11 = factor * s32 * ( 4.0 * dy_yzz - dy_yxx - dy_yyy );
            const double dz_y11 = factor * s32 * ( 4.0 * dz_yzz - dz_yxx - dz_yyy );
            // y12 = factor * ( 2.0 * tgzzz - 3.0 * ( tgzxx + tgzyy ) );
            const double dx_y12 = factor
              * ( 2.0 * dx_zzz - 3.0 * ( dx_zxx + dx_zyy ) );
            const double dy_y12 = factor
              * ( 2.0 * dy_zzz - 3.0 * ( dy_zxx + dy_zyy ) );
            const double dz_y12 = factor
              * ( 2.0 * dz_zzz - 3.0 * ( dz_zxx + dz_zyy ) );
            // y13 = factor * s32 * ( 4.0 * tgxzz - tgxxx - tgxyy );
            const double dx_y13 = factor * s32 * ( 4.0 * dx_xzz - dx_xxx - dx_xyy );
            const double dy_y13 = factor * s32 * ( 4.0 * dy_xzz - dy_xxx - dy_xyy );
            const double dz_y13 = factor * s32 * ( 4.0 * dz_xzz - dz_xxx - dz_xyy );
            // y14 = factor * s15 * ( tgzxx - tgzyy );
            const double dx_y14 = factor * s15 * ( dx_zxx - dx_zyy );
            const double dy_y14 = factor * s15 * ( dy_zxx - dy_zyy );
            const double dz_y14 = factor * s15 * ( dz_zxx - dz_zyy );
            // y15 = factor * s52 * ( tgxxx - 3.0 * tgxyy );
            const double dx_y15 = factor * s52 * ( dx_xxx - 3.0 * dx_xyy );
            const double dy_y15 = factor * s52 * ( dy_xxx - 3.0 * dy_xyy );
            const double dz_y15 = factor * s52 * ( dz_xxx - 3.0 * dz_xyy );

            t09[ig] = y09 * v;
            t10[ig] = y10 * v;
            t11[ig] = y11 * v;
            t12[ig] = y12 * v;
            t13[ig] = y13 * v;
            t14[ig] = y14 * v;
            t15[ig] = y15 * v;

            // contribution to stress tensor
            // sigma_ij = 0.5 * (xi * dj + xj * di)

            dt09_xx[ig] = -( v * dx_y09 * tgx - y09 * dv * tg * tgxx );
            dt09_yy[ig] = -( v * dy_y09 * tgy - y09 * dv * tg * tgyy );
            dt09_zz[ig] = -( v * dz_y09 * tgz - y09 * dv * tg * tgzz );
            dt09_xy[ig] = -( 0.5 * v * ( dx_y09 * tgy + dy_y09 * tgx )
              - y09 * dv * tg * tgxy );
            dt09_yz[ig] = -( 0.5 * v * ( dy_y09 * tgz + dz_y09 * tgy )
              - y09 * dv * tg * tgyz );
            dt09_xz[ig] = -( 0.5 * v * ( dx_y09 * tgz + dz_y09 * tgx )
              - y09 * dv * tg * tgxz );

            dt10_xx[ig] = -( v * dx_y10 * tgx - y10 * dv * tg * tgxx );
            dt10_yy[ig] = -( v * dy_y10 * tgy - y10 * dv * tg * tgyy );
            dt10_zz[ig] = -( v * dz_y10 * tgz - y10 * dv * tg * tgzz );
            dt10_xy[ig] = -( 0.5 * v * ( dx_y10 * tgy + dy_y10 * tgx )
              - y10 * dv * tg * tgxy );
            dt10_yz[ig] = -( 0.5 * v * ( dy_y10 * tgz + dz_y10 * tgy )
              - y10 * dv * tg * tgyz );
            dt10_xz[ig] = -( 0.5 * v * ( dx_y10 * tgz + dz_y10 * tgx )
              - y10 * dv * tg * tgxz );

            dt11_xx[ig] = -( v * dx_y11 * tgx - y11 * dv * tg * tgxx );
            dt11_yy[ig] = -( v * dy_y11 * tgy - y11 * dv * tg * tgyy );
            dt11_zz[ig] = -( v * dz_y11 * tgz - y11 * dv * tg * tgzz );
            dt11_xy[ig] = -( 0.5 * v * ( dx_y11 * tgy + dy_y11 * tgx )
              - y11 * dv * tg * tgxy );
            dt11_yz[ig] = -( 0.5 * v * ( dy_y11 * tgz + dz_y11 * tgy )
              - y11 * dv * tg * tgyz );
            dt11_xz[ig] = -( 0.5 * v * ( dx_y11 * tgz + dz_y11 * tgx )
              - y11 * dv * tg * tgxz );

            dt12_xx[ig] = -( v * dx_y12 * tgx - y12 * dv * tg * tgxx );
            dt12_yy[ig] = -( v * dy_y12 * tgy - y12 * dv * tg * tgyy );
            dt12_zz[ig] = -( v * dz_y12 * tgz - y12 * dv * tg * tgzz );
            dt12_xy[ig] = -( 0.5 * v * ( dx_y12 * tgy + dy_y12 * tgx )
              - y12 * dv * tg * tgxy );
            dt12_yz[ig] = -( 0.5 * v * ( dy_y12 * tgz + dz_y12 * tgy )
              - y12 * dv * tg * tgyz );
            dt12_xz[ig] = -( 0.5 * v * ( dx_y12 * tgz + dz_y12 * tgx )
              - y12 * dv * tg * tgxz );

            dt13_xx[ig] = -( v * dx_y13 * tgx - y13 * dv * tg * tgxx );
            dt13_yy[ig] = -( v * dy_y13 * tgy - y13 * dv * tg * tgyy );
            dt13_zz[ig] = -( v * dz_y13 * tgz - y13 * dv * tg * tgzz );
            dt13_xy[ig] = -( 0.5 * v * ( dx_y13 * tgy + dy_y13 * tgx )
              - y13 * dv * tg * tgxy );
            dt13_yz[ig] = -( 0.5 * v * ( dy_y13 * tgz + dz_y13 * tgy )
              - y13 * dv * tg * tgyz );
            dt13_xz[ig] = -( 0.5 * v * ( dx_y13 * tgz + dz_y13 * tgx )
              - y13 * dv * tg * tgxz );

            dt14_xx[ig] = -( v * dx_y14 * tgx - y14 * dv * tg * tgxx );
            dt14_yy[ig] = -( v * dy_y14 * tgy - y14 * dv * tg * tgyy );
            dt14_zz[ig] = -( v * dz_y14 * tgz - y14 * dv * tg * tgzz );
            dt14_xy[ig] = -( 0.5 * v * ( dx_y14 * tgy + dy_y14 * tgx )
              - y14 * dv * tg * tgxy );
            dt14_yz[ig] = -( 0.5 * v * ( dy_y14 * tgz + dz_y14 * tgy )
              - y14 * dv * tg * tgyz );
            dt14_xz[ig] = -( 0.5 * v * ( dx_y14 * tgz + dz_y14 * tgx )
              - y14 * dv * tg * tgxz );

            dt15_xx[ig] = -( v * dx_y15 * tgx - y15 * dv * tg * tgxx );
            dt15_yy[ig] = -( v * dy_y15 * tgy - y15 * dv * tg * tgyy );
            dt15_zz[ig] = -( v * dz_y15 * tgz - y15 * dv * tg * tgzz );
            dt15_xy[ig] = -( 0.5 * v * ( dx_y15 * tgy + dy_y15 * tgx )
              - y15 * dv * tg * tgxy );
            dt15_yz[ig] = -( 0.5 * v * ( dy_y15 * tgz + dz_y15 * tgy )
              - y15 * dv * tg * tgyz );
            dt15_xz[ig] = -( 0.5 * v * ( dx_y15 * tgz + dz_y15 * tgx )
              - y15 * dv * tg * tgxz );

          }

          ilm += 2 * l + 1;

        } // l == 3
        else
        {
          assert(false);
        }
      } // l != lloc[is]
    } // for l
    assert(ilm == s->nlm());
  }
  tmap["update_twnl"].stop();
}

////////////////////////////////////////////////////////////////////////////////
double NonLocalPotential::energy(bool compute_hpsi, SlaterDet& dsd,
    bool compute_forces, vector<vector<double> >& fion_enl,
    bool compute_stress, valarray<double>& sigma_enl)
{
  for ( int is = 0; is < fion_enl.size(); is++ )
    for ( int i = 0; i < fion_enl[is].size(); i++ )
       fion_enl[is][i] = 0.0;
  sigma_enl = 0.0;

  if ( nspnl == 0 ) return 0.0;

  double enl = 0.0;
  double tsum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const vector<double>& occ = sd_.occ();
  const int ngwl = basis_.localsize();
  // define atom block size
  const int na_block_size = 32;
  valarray<double> gr(na_block_size*ngwl); // gr[ig+ia*ngwl]
  valarray<double> cgr(na_block_size*ngwl); // cgr[ig+ia*ngwl]
  valarray<double> sgr(na_block_size*ngwl); // sgr[ig+ia*ngwl]
  vector<vector<double> > tau;
  atoms_.get_positions(tau);

  const double omega = basis_.cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  for ( int is = 0; is < nsp; is++ )
  {
    Species *s = atoms_.species_list[is];
    if ( npr[is] > 0 ) // species is is non-local
    {
      valarray<double> tmpfion(3*na[is]);
      tmpfion = 0.0;
      // define number of atom blocks
      const int na_blocks = na[is] / na_block_size +
                            ( na[is] % na_block_size == 0 ? 0 : 1 );

      valarray<double> anl_loc(npr[is]*na_block_size*2*ngwl);
      const int nstloc = sd_.nstloc();
      // fnl_loc[ipra][n]
      // fnl is real if basis is real, complex otherwise
      const int fnl_loc_size = basis_.real() ? npr[is]*na_block_size*nstloc :
        2*npr[is]*na_block_size*nstloc;
      valarray<double> fnl_loc(fnl_loc_size);
      valarray<double> fnl_buf(fnl_loc_size);
      for ( int ia_block = 0; ia_block < na_blocks; ia_block++ )
      {
        // process projectors of atoms in block ia_block

        const int iastart = ia_block * na_block_size;
        const int iaend = (ia_block+1) * na_block_size < na[is] ?
                        (ia_block+1) * na_block_size :
                        na[is];
        const int ia_block_size = iaend - iastart;

        // compute cgr[is][ia][ig], sgr[is][ia][ig]
        tmap["comp_eigr"].start();
        int k = 3;
        double mone = -1.0, zero = 0.0;
        char cn='n';

        // next line: const cast is ok since dgemm_ does not modify argument
        double* kpgx = const_cast<double*>(basis_.kpgx_ptr(0));
        dgemm(&cn,&cn,(int*)&ngwl,(int*)&ia_block_size,&k,&mone,
              kpgx,(int*)&ngwl, &tau[is][3*iastart],&k,
              &zero,&gr[0],(int*)&ngwl);

        int len = ia_block_size * ngwl;
#if USE_MASSV
        vsincos(&sgr[0],&cgr[0],&gr[0],&len);
#else
        for ( int i = 0; i < len; i++ )
        {
          const double arg = gr[i];
          sgr[i] = sin(arg);
          cgr[i] = cos(arg);
        }
#endif
        tmap["comp_eigr"].stop();

        // compute anl_loc
        tmap["comp_anl"].start();
        for ( int ipr = 0; ipr < npr[is]; ipr++ )
        {
          // twnl[is][ig+ngwl*ipr]
          const double * t = &twnl[is][ngwl*ipr];
          const int l = lproj[is][ipr];

          // anl_loc[ig+ipra*ngwl]

          for ( int ia = 0; ia < ia_block_size; ia++ )
          {
            double* a = &anl_loc[2*(ia+ipr*ia_block_size)*ngwl];
            const double* c = &cgr[ia*ngwl];
            const double* s = &sgr[ia*ngwl];
            if ( l == 0 )
            {
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                a[2*ig]   = t[ig] * c[ig];
                a[2*ig+1] = t[ig] * s[ig];
              }
            }
            else if ( l == 1 )
            {
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                /* Next line: -i * eigr */
                /* -i * (a+i*b) = b - i*a */
                a[2*ig]   =  t[ig] * s[ig];
                a[2*ig+1] = -t[ig] * c[ig];
              }
            }
            else if ( l == 2 )
            {
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                // Next line: (-) sign for -eigr
                a[2*ig]   = -t[ig] * c[ig];
                a[2*ig+1] = -t[ig] * s[ig];
              }
            }
            else if ( l == 3 )
            {
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                // Next line: i * eigr
                // i * (a+i*b) = -b + i*a
                a[2*ig]   = -t[ig] * s[ig];
                a[2*ig+1] =  t[ig] * c[ig];
              }
            }
          }
        } // ipr
        tmap["comp_anl"].stop();

        // array anl_loc is complete

        // compute fnl[npra][nstloc] = anl^T * c
        int nprnaloc = ia_block_size * npr[is];
        if ( basis_.real() )
        {
          double one=1.0;
          char ct='t';
          int twongwl = 2 * ngwl;
          int c_lda = 2*sd_.c().mloc();
          const complex<double>* c = sd_.c().cvalptr();
          tmap["fnl_gemm"].start();
          dgemm(&ct,&cn,&nprnaloc,(int*)&nstloc,&twongwl,&one,
                &anl_loc[0],&twongwl, (double*)c, &c_lda,
                &zero,&fnl_loc[0],&nprnaloc);
          tmap["fnl_gemm"].stop();

          // correct for double counting if ctxt_.myrow() == 0
          if ( ctxt_.myrow() == 0 )
          {
            // rank-one update
            // dger(m,n,alpha,x,incx,y,incy,a,lda);
            // a += alpha * x * transpose(y)
            // x = first row of anl_loc
            // y^T = first row of c
            double alpha = -0.5;
            dger(&nprnaloc,(int*)&nstloc,&alpha,&anl_loc[0],&twongwl,
                 (double*)c,&c_lda,&fnl_loc[0],&nprnaloc);
          }
        }
        else
        {
          // fnl is complex
          complex<double> cone=1.0;
          complex<double> czero=0.0;
          char cc='c';
          int c_lda = sd_.c().mloc();
          const complex<double>* c = sd_.c().cvalptr();
          tmap["fnl_zemm"].start();
          zgemm(&cc,&cn,&nprnaloc,(int*)&nstloc,(int*)&ngwl,&cone,
                (complex<double>*) &anl_loc[0],(int*)&ngwl,
                (complex<double>*) c, &c_lda,&czero,
                (complex<double>*)&fnl_loc[0],&nprnaloc);
          tmap["fnl_zemm"].stop();
        }

#if USE_MPI
        tmap["fnl_allreduce"].start();
        // Allreduce fnl partial sum
        MPI_Comm basis_comm = basis_.comm();
        int fnl_size = basis_.real() ? nprnaloc*nstloc : 2*nprnaloc*nstloc;
        MPI_Allreduce(&fnl_loc[0],&fnl_buf[0],fnl_size,
                      MPI_DOUBLE,MPI_SUM,basis_comm);
        tmap["fnl_allreduce"].stop();

        // factor 2.0 in next line is: counting G, -G
        if ( basis_.real() )
          fnl_loc = 2.0 * fnl_buf;
        else
          fnl_loc = fnl_buf;
#else
        // factor 2.0 in next line is: counting G, -G
        if ( basis_.real() )
          fnl_loc *= 2.0;
#endif

        // if the species has multiple projectors, that are not orthogonal
        // multiply with D matrix
        if ( s->has_dmatrix() )
        {
          if ( basis_.real() )
          {
            // helper variables
            double one = 1.0;
            double zero = 0.0;
            const size_t dsize = nprnaloc * nprnaloc;
            //
            // construct D matrix
            //
            // allocate array
            valarray<double> dmatrix(zero,dsize);
            // loop over all elements
            for ( int i = 0; i < nprnaloc; i++ )
            {
              // extract atom, l, and m corresponding to this index
              const int iatom = i % ia_block_size;
              const int ipr = i / ia_block_size;
              for ( int j = 0; j <= i; j++ )
              {
                // extract atom, l, and m corresponding to this index
                const int jatom = j % ia_block_size;
                const int jpr = j / ia_block_size;
                // D matrix is diagonal in atom, l, and m
                if ( iatom != jatom ) continue;
                const double val = s->dmatrix(ipr,jpr);
                dmatrix[nprnaloc * i + j] = val;
                // matrix is symmetric
                if ( i != j ) dmatrix[nprnaloc * j + i] = val;
              }
            }
            // multiply F'_In = D_IJ F_Jn
            dgemm(&cn,&cn,&nprnaloc,(int*) &nstloc,&nprnaloc,&one,&dmatrix[0],
              &nprnaloc,&fnl_loc[0],&nprnaloc,&zero,&fnl_buf[0],&nprnaloc);
          }
          else
          {
            // helper variables
            complex<double> one = 1.0;
            complex<double> zero = 0.0;
            const size_t dsize = nprnaloc * nprnaloc;
            //
            // construct D matrix
            //
            // allocate array
            valarray < complex<double> > dmatrix(zero,dsize);
            // loop over all elements
            for ( int i = 0; i < nprnaloc; i++ )
            {
              // extract atom, l, and m corresponding to this index
              const int iatom = i % ia_block_size;
              const int ipr = i / ia_block_size;
              for ( int j = 0; j <= i; j++ )
              {
                // extract atom, l, and m corresponding to this index
                const int jatom = j % ia_block_size;
                const int jpr = j / ia_block_size;
                // D matrix is diagonal in atom, l, and m
                if ( iatom != jatom ) continue;
                const double val = s->dmatrix(ipr,jpr);
                dmatrix[nprnaloc * i + j] = val;
                // matrix is symmetric
                if ( i != j ) dmatrix[nprnaloc * j + i] = val;
              }
            }
            // multiply F'_In = D_IJ F_Jn
            zgemm(&cn,&cn,&nprnaloc,(int*) &nstloc,&nprnaloc,&one,&dmatrix[0],
              &nprnaloc,(complex<double>*) &fnl_loc[0],&nprnaloc,&zero,
              (complex<double>*) &fnl_buf[0],&nprnaloc);
          }
        }
        else
        {
          fnl_buf = fnl_loc;
        }

        // accumulate Enl contribution
        const int nbase = ctxt_.mycol() * sd_.c().nb();
        if ( basis_.real() )
        {
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            const double fac = wt[is][ipr] * omega_inv;
            for ( int n = 0; n < nstloc; n++ )
            {
              const double facn = fac * occ[n + nbase];
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                const int i = ia + ipr*ia_block_size + n * nprnaloc;
                //cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia
                //     << ",n=" << n << "]: " << fnl_loc[i] << endl;
                enl += facn * fnl_loc[i] * fnl_buf[i];
                fnl_loc[i] = fac * fnl_buf[i];
              }
            }
          }
        }
        else
        {
          // fnl is complex
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            const double fac = wt[is][ipr] * omega_inv;
            for ( int n = 0; n < nstloc; n++ )
            {
              const double facn = fac * occ[n + nbase];
              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                const int i = ia + ipr*ia_block_size + n * nprnaloc;
                //cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia
                //     << ",n=" << n << "]: " << fnl_loc[i] << endl;
                const double f_re = fnl_loc[2*i];
                const double f_im = fnl_loc[2*i+1];
                const double fb_re = fnl_buf[2*i];
                const double fb_im = fnl_buf[2*i+1];
                enl += facn * (f_re*fb_re + f_im*fb_im);
                fnl_loc[2*i] = fac * fb_re;
                fnl_loc[2*i+1] = fac * fb_im;
              }
            }
          }
        }

        if ( compute_hpsi )
        {
          if ( basis_.real() )
          {
            tmap["enl_hpsi"].start();
            // compute cp += anl * fnl
            complex<double>* cp = dsd.c().valptr();
            int cp_lda = 2*dsd.c().mloc();
            int twongwl = 2 * ngwl;
            double one = 1.0;
            dgemm(&cn,&cn,&twongwl,(int*)&nstloc,&nprnaloc,&one,
                  &anl_loc[0],&twongwl, &fnl_loc[0],&nprnaloc,
                  &one,(double*)cp, &cp_lda);
            tmap["enl_hpsi"].stop();
          }
          else
          {
            tmap["enl_hpsi"].start();
            // compute cp += anl * fnl
            complex<double>* cp = dsd.c().valptr();
            complex<double> cone = 1.0;
            int cp_lda = dsd.c().mloc();
            zgemm(&cn,&cn,(int*)&ngwl,(int*)&nstloc,&nprnaloc,&cone,
                  (complex<double>*) &anl_loc[0],(int*)&ngwl,
                  (complex<double>*) &fnl_loc[0],&nprnaloc,
                  &cone,cp, &cp_lda);
            tmap["enl_hpsi"].stop();
          }
        }

        // ionic forces
        if ( compute_forces )
        {
          tmap["enl_fion"].start();

          valarray<double> dfnl_loc(fnl_loc_size);
          for ( int j = 0; j < 3; j++ )
          {
            const double *const kpgxj = basis_.kpgx_ptr(j);

            // compute anl_loc
            for ( int ipr = 0; ipr < npr[is]; ipr++ )
            {
              // twnl[is][ig+ngwl*ipr]
              const double * t = &twnl[is][ngwl*ipr];
              const int l = lproj[is][ipr];

              // anl_loc[ig+ipra*ngwl]

              for ( int ia = 0; ia < ia_block_size; ia++ )
              {
                double* a = &anl_loc[2*(ia+ipr*ia_block_size)*ngwl];
                const double* c = &cgr[ia*ngwl];
                const double* s = &sgr[ia*ngwl];
                if ( l == 0 )
                {
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    const double tt = kpgxj[ig] * t[ig];
                    // Next lines: -i * ( a + ib ) = b - ia
                    a[2*ig]   =  tt * s[ig];
                    a[2*ig+1] = -tt * c[ig];
                  }
                }
                else if ( l == 1 )
                {
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    // Next lines: (-i)**2 * ( a + ib ) = - a - ib
                    const double tt = - kpgxj[ig] * t[ig];
                    a[2*ig]   = tt * c[ig];
                    a[2*ig+1] = tt * s[ig];
                  }
                }
                else if ( l == 2 )
                {
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    // Next lines: (-i) * - ( a + ib ) = i*(a+ib) = - b + ia
                    const double tt = kpgxj[ig] * t[ig];
                    a[2*ig]   = -tt * s[ig];
                    a[2*ig+1] =  tt * c[ig];
                  }
                }
                else if ( l == 3 )
                {
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    // Next lines: (-i)**4 * ( a + ib ) = a + ib
                    const double tt = kpgxj[ig] * t[ig];
                    a[2*ig]   = tt * c[ig];
                    a[2*ig+1] = tt * s[ig];
                  }
                }
              }
            } // ipr

            // array anl_loc is complete

            // compute dfnl[npra][nstloc] = anl^T * c
            if ( basis_.real() )
            {
              // factor 2.0 in next line is: counting G, -G
              // Note: no need to correct for double counting of the
              // G=0 component which is always zero
              double two=2.0;
              char ct='t';
              const int twongwl = 2 * ngwl;
              const int nprnaloc = ia_block_size * npr[is];
              int c_lda = 2*sd_.c().mloc();
              const complex<double>* c = sd_.c().cvalptr();
              dgemm(&ct,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&twongwl,&two,
                    &anl_loc[0],(int*)&twongwl, (double*)c,(int*)&c_lda,
                    &zero,&dfnl_loc[0],(int*)&nprnaloc);
            }
            else
            {
              complex<double> cone=1.0;
              complex<double> czero=0.0;
              char cc='c';
              const int nprnaloc = ia_block_size * npr[is];
              int c_lda = sd_.c().mloc();
              const complex<double>* c = sd_.c().cvalptr();
              zgemm(&cc,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&ngwl,&cone,
                    (complex<double>*) &anl_loc[0],(int*)&ngwl,
                    (complex<double>*) c,(int*)&c_lda,
                    &czero,(complex<double>*)&dfnl_loc[0],(int*)&nprnaloc);
            }

            // accumulate non-local contributions to forces
            if ( basis_.real() )
            {
              for ( int ipr = 0; ipr < npr[is]; ipr++ )
              {
                for ( int n = 0; n < nstloc; n++ )
                {
                  // Factor 2.0 in next line from derivative of |Fnl|^2
                  const double facn = 2.0 * occ[n + nbase];
                  for ( int ia = 0; ia < ia_block_size; ia++ )
                  {
                    const int ia_global = ia + iastart;
                    const int i = ia + ipr*ia_block_size + n * nprnaloc;
                    //cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia
                    //     << ",n=" << n << "]: " << fnl_loc[i] << endl;
                    tmpfion[3*ia_global+j] -= facn *
                                              fnl_loc[i] * dfnl_loc[i];
                  }
                }
              }
            }
            else
            {
              for ( int ipr = 0; ipr < npr[is]; ipr++ )
              {
                for ( int n = 0; n < nstloc; n++ )
                {
                  // Factor 2.0 in next line from derivative of |Fnl|^2
                  const double facn = 2.0 * occ[n + nbase];
                  for ( int ia = 0; ia < ia_block_size; ia++ )
                  {
                    const int ia_global = ia + iastart;
                    const int i = ia + ipr*ia_block_size + n * nprnaloc;
                    //cout << "fnl_loc[ipr=" << ipr << ",ia=" << ia
                    //     << ",n=" << n << "]: " << fnl_loc[i] << endl;
                    double f_re = fnl_loc[2*i];
                    double f_im = fnl_loc[2*i+1];
                    double df_re = dfnl_loc[2*i];
                    double df_im = dfnl_loc[2*i+1];
                    tmpfion[3*ia_global+j] -=
                      facn * ( f_re * df_re + f_im * df_im );
                  }
                }
              }
            }
          } // j

          tmap["enl_fion"].stop();
        } // compute_forces

        if ( compute_stress )
        {
          tmap["enl_sigma"].start();
          valarray<double> dfnl_loc(fnl_loc_size);

          for ( int ij = 0; ij < 6; ij++ )
          {
            // compute anl_loc
            int ipr = 0;
            while ( ipr < npr[is] )
            {
              // twnl[is][ig+ngwl*ipr]
              const int l = lproj[is][ipr];
              if ( l == 0 )
              {
                // dtwnl[is][ipr][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * ipr))
                // ipr = iquad + nquad[is] * ilm, where ilm = 0
                const double *const dt0 = &dtwnl[is][ngwl*(ij+6*ipr)];
                for ( int ia = 0; ia < ia_block_size; ia++ )
                {
                  double* a0 = &anl_loc[2*(ia+ipr*ia_block_size)*ngwl];
                  const double* c = &cgr[ia*ngwl];
                  const double* s = &sgr[ia*ngwl];
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    const double d0 = dt0[ig];
                    // Next lines: -i * ( a + ib ) = b - ia
                    a0[2*ig]   =  d0 * c[ig];
                    a0[2*ig+1] =  d0 * s[ig];
                  }
                }
              }
              else if ( l == 1 )
              {
                const int ipr1 = ipr;
                const int ipr2 = ipr + 1;
                const int ipr3 = ipr + 2;
                // dtwnl[is][ipr][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * iprx ))
                const double *dt1 = &dtwnl[is][ngwl*(ij+6*ipr1)];
                const double *dt2 = &dtwnl[is][ngwl*(ij+6*ipr2)];
                const double *dt3 = &dtwnl[is][ngwl*(ij+6*ipr3)];
                for ( int ia = 0; ia < ia_block_size; ia++ )
                {
                  double* a1 = &anl_loc[2*(ia+ipr1*ia_block_size)*ngwl];
                  double* a2 = &anl_loc[2*(ia+ipr2*ia_block_size)*ngwl];
                  double* a3 = &anl_loc[2*(ia+ipr3*ia_block_size)*ngwl];
                  const double* c = &cgr[ia*ngwl];
                  const double* s = &sgr[ia*ngwl];
                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    const double d1 = dt1[ig];
                    const double d2 = dt2[ig];
                    const double d3 = dt3[ig];
                    // Next line: (-i)^l factor is -i
                    // Next line: -i * eigr
                    // -i * (a+i*b) = b - i*a
                    const double tc = -c[ig]; //  -cosgr[ia][ig]
                    const double ts =  s[ig]; //   singr[ia][ig]
                    a1[2*ig]   = d1 * ts;
                    a1[2*ig+1] = d1 * tc;
                    a2[2*ig]   = d2 * ts;
                    a2[2*ig+1] = d2 * tc;
                    a3[2*ig]   = d3 * ts;
                    a3[2*ig+1] = d3 * tc;
                  }
                }
              }
              else if ( l == 2 )
              {
                const int ipr4 = ipr;
                const int ipr5 = ipr + 1;
                const int ipr6 = ipr + 2;
                const int ipr7 = ipr + 3;
                const int ipr8 = ipr + 4;
                // dtwnl[is][ipr][iquad][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * ( iquad + nquad[is] * ipr ))
                const double *dt4 = &dtwnl[is][ngwl*(ij+6*ipr4)];
                const double *dt5 = &dtwnl[is][ngwl*(ij+6*ipr5)];
                const double *dt6 = &dtwnl[is][ngwl*(ij+6*ipr6)];
                const double *dt7 = &dtwnl[is][ngwl*(ij+6*ipr7)];
                const double *dt8 = &dtwnl[is][ngwl*(ij+6*ipr8)];
                for ( int ia = 0; ia < ia_block_size; ia++ )
                {
                  double* a4 = &anl_loc[2*(ia+ipr4*ia_block_size)*ngwl];
                  double* a5 = &anl_loc[2*(ia+ipr5*ia_block_size)*ngwl];
                  double* a6 = &anl_loc[2*(ia+ipr6*ia_block_size)*ngwl];
                  double* a7 = &anl_loc[2*(ia+ipr7*ia_block_size)*ngwl];
                  double* a8 = &anl_loc[2*(ia+ipr8*ia_block_size)*ngwl];
                  const double* c = &cgr[ia*ngwl];
                  const double* s = &sgr[ia*ngwl];

                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    const double d4 = dt4[ig];
                    const double d5 = dt5[ig];
                    const double d6 = dt6[ig];
                    const double d7 = dt7[ig];
                    const double d8 = dt8[ig];
                    // Next lines: (-i)^2 * ( a + ib ) =  - ( a + ib )
                    const double tc = -c[ig]; //  -cosgr[ia][ig]
                    const double ts = -s[ig]; //  -singr[ia][ig]
                    a4[2*ig]   = d4 * tc;
                    a4[2*ig+1] = d4 * ts;
                    a5[2*ig]   = d5 * tc;
                    a5[2*ig+1] = d5 * ts;
                    a6[2*ig]   = d6 * tc;
                    a6[2*ig+1] = d6 * ts;
                    a7[2*ig]   = d7 * tc;
                    a7[2*ig+1] = d7 * ts;
                    a8[2*ig]   = d8 * tc;
                    a8[2*ig+1] = d8 * ts;
                  }
                }
              }
              else if ( l == 3 )
              {
                const int ipr09 = ipr;
                const int ipr10= ipr + 1;
                const int ipr11= ipr + 2;
                const int ipr12= ipr + 3;
                const int ipr13= ipr + 4;
                const int ipr14= ipr + 5;
                const int ipr15= ipr + 6;
                // dtwnl[is][ipr][iquad][ij][ngwl]
                // index = ig + ngwl * ( ij + 6 * ( iquad + nquad[is] * ipr ))
                const double *dt09 = &dtwnl[is][ngwl * ( ij + 6 * ipr09 )];
                const double *dt10 = &dtwnl[is][ngwl * ( ij + 6 * ipr10 )];
                const double *dt11 = &dtwnl[is][ngwl * ( ij + 6 * ipr11 )];
                const double *dt12 = &dtwnl[is][ngwl * ( ij + 6 * ipr12 )];
                const double *dt13 = &dtwnl[is][ngwl * ( ij + 6 * ipr13 )];
                const double *dt14 = &dtwnl[is][ngwl * ( ij + 6 * ipr14 )];
                const double *dt15 = &dtwnl[is][ngwl * ( ij + 6 * ipr15 )];
                for ( int ia = 0; ia < ia_block_size; ia++ )
                {
                  double* a09 = &anl_loc[2 * ( ia + ipr09 * ia_block_size ) * ngwl];
                  double* a10 = &anl_loc[2 * ( ia + ipr10 * ia_block_size ) * ngwl];
                  double* a11 = &anl_loc[2 * ( ia + ipr11 * ia_block_size ) * ngwl];
                  double* a12 = &anl_loc[2 * ( ia + ipr12 * ia_block_size ) * ngwl];
                  double* a13 = &anl_loc[2 * ( ia + ipr13 * ia_block_size ) * ngwl];
                  double* a14 = &anl_loc[2 * ( ia + ipr14 * ia_block_size ) * ngwl];
                  double* a15 = &anl_loc[2 * ( ia + ipr15 * ia_block_size ) * ngwl];
                  const double* c = &cgr[ia * ngwl];
                  const double* s = &sgr[ia * ngwl];

                  for ( int ig = 0; ig < ngwl; ig++ )
                  {
                    const double d09 = dt09[ig];
                    const double d10 = dt10[ig];
                    const double d11 = dt11[ig];
                    const double d12 = dt12[ig];
                    const double d13 = dt13[ig];
                    const double d14 = dt14[ig];
                    const double d15 = dt15[ig];
                    // Next lines: (-i)^2 * ( a + ib ) =  - ( a + ib )
                    const double tc =  c[ig]; //   cosgr[ia][ig]
                    const double ts = -s[ig]; //  -singr[ia][ig]
                    a09[2 * ig] = d09 * ts;
                    a09[2 * ig + 1] = d09 * tc;
                    a10[2 * ig] = d10 * ts;
                    a10[2 * ig + 1] = d10 * tc;
                    a11[2 * ig] = d11 * ts;
                    a11[2 * ig + 1] = d11 * tc;
                    a12[2 * ig] = d12 * ts;
                    a12[2 * ig + 1] = d12 * tc;
                    a13[2 * ig] = d13 * ts;
                    a13[2 * ig + 1] = d13 * tc;
                    a14[2 * ig] = d14 * ts;
                    a14[2 * ig + 1] = d14 * tc;
                    a15[2 * ig] = d15 * ts;
                    a15[2 * ig + 1] = d15 * tc;
                  }
                }
              }
              else
              {
                assert(false);
              } // l
              ipr += 2*l+1;
            } // while ipr

            // array anl_loc is complete

            // compute dfnl[npra][nstloc] = anl^T * c
            if ( basis_.real() )
            {
              // factor 2.0 in next line is: counting G, -G
              // Note: no need to correct for double counting of the
              // G=0 component which is always zero
              double two=2.0;
              char ct='t';
              const int twongwl = 2 * ngwl;
              const int nprnaloc = ia_block_size * npr[is];
              int c_lda = 2*sd_.c().mloc();
              const complex<double>* c = sd_.c().cvalptr();
              dgemm(&ct,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&twongwl,&two,
                    &anl_loc[0],(int*)&twongwl, (double*)c,(int*)&c_lda,
                    &zero,&dfnl_loc[0],(int*)&nprnaloc);
            }
            else
            {
              complex<double> cone=1.0;
              complex<double> czero=0.0;
              char cc='c';
              const int nprnaloc = ia_block_size * npr[is];
              int c_lda = sd_.c().mloc();
              const complex<double>* c = sd_.c().cvalptr();
              zgemm(&cc,&cn,(int*)&nprnaloc,(int*)&nstloc,(int*)&ngwl,&cone,
                    (complex<double>*)&anl_loc[0],(int*)&ngwl,
                    (complex<double>*)c,(int*)&c_lda,
                    &czero,(complex<double>*)&dfnl_loc[0],(int*)&nprnaloc);
            }

            // accumulate non-local contributions to sigma_ij
            if ( basis_.real() )
            {
              for ( int n = 0; n < nstloc; n++ )
              {
                // Factor 2.0 in next line from derivative of |Fnl|^2
                const double facn = 2.0 * occ[n + nbase];
                for ( int ipra = 0; ipra < npr[is]*ia_block_size; ipra++ )
                {
                  const int i = ipra + n * nprnaloc;
                  tsum[ij] += facn * fnl_loc[i] * dfnl_loc[i];
                }
              }
            }
            else
            {
              for ( int n = 0; n < nstloc; n++ )
              {
                // Factor 2.0 in next line from derivative of |Fnl|^2
                const double facn = 2.0 * occ[n + nbase];
                for ( int ipra = 0; ipra < npr[is]*ia_block_size; ipra++ )
                {
                  const int i = ipra + n * nprnaloc;
                  double f_re = fnl_loc[2*i];
                  double f_im = fnl_loc[2*i+1];
                  double df_re = dfnl_loc[2*i];
                  double df_im = dfnl_loc[2*i+1];
                  tsum[ij] += facn * ( f_re * df_re + f_im * df_im );
                }
              }
            }
          } // ij
          tmap["enl_sigma"].stop();
        } // compute_stress

      } // ia_block

      if ( compute_forces )
      {
        ctxt_.dsum(3*na[is],1,&tmpfion[0],3*na[is]);
        for ( int ia = 0; ia < na[is]; ia++ )
        {
          fion_enl[is][3*ia+0] += tmpfion[3*ia];
          fion_enl[is][3*ia+1] += tmpfion[3*ia+1];
          fion_enl[is][3*ia+2] += tmpfion[3*ia+2];
        }
      }
    } // npr[is]>0
  } // is

  // reduction of enl across rows
  ctxt_.dsum('r',1,1,&enl,1);

  if ( compute_stress )
  {
    ctxt_.dsum(6,1,&tsum[0],6);
    sigma_enl[0] = ( enl + tsum[0] ) * omega_inv;
    sigma_enl[1] = ( enl + tsum[1] ) * omega_inv;
    sigma_enl[2] = ( enl + tsum[2] ) * omega_inv;
    sigma_enl[3] = + tsum[3] * omega_inv;
    sigma_enl[4] = + tsum[4] * omega_inv;
    sigma_enl[5] = + tsum[5] * omega_inv;
  }

  return enl;
}
