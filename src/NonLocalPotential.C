////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.C,v 1.23 2007-10-19 17:15:23 fgygi Exp $

#include "NonLocalPotential.h"
#include "Species.h"
#include "blas.h"
#include <iomanip>
using namespace std;

#if AIX || BGL
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
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::init(void)
{
  const int ngwl = basis_.localsize();

  nsp = atoms_.nsp();

  lmax.resize(nsp);
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
      lmax[is] = s->lmax();
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
      for ( int l = 0; l <= lmax[is]; l++ )
      {
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
              wt[is][ipr] = s->wsg(l);
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

  const int ngwl = basis_.localsize();
  const double pi = M_PI;
  const double fpi = 4.0 * pi;
  const double s14pi = sqrt(1.0/fpi);
  const double s34pi = sqrt(3.0/fpi);
  const double s54pi = sqrt(5.0/fpi);
  const double s20pi = sqrt(20.0*pi);
  const double s20pi3 = sqrt(20.0*pi/3.0);
  const double s3 = sqrt(3.0);

  const double *kpg   = basis_.kpg_ptr();
  const double *kpg2  = basis_.kpg2_ptr();
  const double *kpgi  = basis_.kpgi_ptr();
  const double *kpg_x = basis_.kpgx_ptr(0);
  const double *kpg_y = basis_.kpgx_ptr(1);
  const double *kpg_z = basis_.kpgx_ptr(2);

  // compute twnl and dtwnl
  for ( int is = 0; is < nsp; is++ )
  {
    Species *s = atoms_.species_list[is];

    int ilm = 0;
    for ( int l = 0; l <= lmax[is]; l++ )
    {
      if ( l != lloc[is] )
      {
        if ( l == 0 )
        {
          if ( nquad[is] == 0 )
          {
            // Kleinman-Bylander

            // twnl[is][ipr][ig]
            // ipr = ilm = 0
            // index = ig + ngwl*ipr, i.e. index = ig
            double *t0  = &twnl[is][0];

            // dtwnl[is][ipr][ij][ngwl]
            // index = ig + ngwl * ( ij + 6 * ipr ), ipr = 0
            // i.e. index = ig + ij * ngwl
            double *dt0_xx = &dtwnl[is][0*ngwl];
            double *dt0_yy = &dtwnl[is][1*ngwl];
            double *dt0_zz = &dtwnl[is][2*ngwl];
            double *dt0_xy = &dtwnl[is][3*ngwl];
            double *dt0_yz = &dtwnl[is][4*ngwl];
            double *dt0_xz = &dtwnl[is][5*ngwl];
            // Special case k=G=0 is ok since kpgi[0] = 0.0 at k=G=0
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              s->dvnlg(0,kpg[ig],v,dv);

              t0[ig] = s14pi * v;

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;

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
              s->dvnlg(l,tg,v,dv);

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

              s->dvnlg(l,tg,v,dv);

              const double tgx = kpg_x[ig];
              const double tgy = kpg_y[ig];
              const double tgz = kpg_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;

              const double tgi = kpgi[ig];
              const double tg2 = tg * tg;
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
                const double tg2 = tg * tg;
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
    bool compute_forces, vector<vector<double> >& fion,
    bool compute_stress, valarray<double>& sigma_enl)
{
  const vector<double>& occ = sd_.occ();
  const int ngwl = basis_.localsize();
  // define atom block size
  const int na_block_size = 32;
  valarray<double> gr(na_block_size*ngwl); // gr[ig+ia*ngwl]
  valarray<double> cgr(na_block_size*ngwl); // cgr[ig+ia*ngwl]
  valarray<double> sgr(na_block_size*ngwl); // sgr[ig+ia*ngwl]
  vector<vector<double> > tau;
  atoms_.get_positions(tau);

  double enl = 0.0;
  double tsum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  if ( nspnl == 0 ) return 0.0;
  const double omega = basis_.cell().volume();
  assert(omega != 0.0);
  const double omega_inv = 1.0 / omega;

  for ( int is = 0; is < nsp; is++ )
  {
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
#if AIX || BGL
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

        tmap["fnl_allreduce"].start();
        // Allreduce fnl partial sum
        MPI_Comm basis_comm = basis_.context().comm();
        int fnl_size = basis_.real() ? nprnaloc*nstloc : 2*nprnaloc*nstloc;
        MPI_Allreduce(&fnl_loc[0],&fnl_buf[0],fnl_size,
                      MPI_DOUBLE,MPI_SUM,basis_comm);
        tmap["fnl_allreduce"].stop();

        // factor 2.0 in next line is: counting G, -G
        if ( basis_.real() )
          fnl_loc = 2.0 * fnl_buf;
        else
          fnl_loc = fnl_buf;

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
                const double tmp = fnl_loc[i];
                enl += facn * tmp * tmp;
                fnl_loc[i] = fac * tmp;
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
                enl += facn * (f_re*f_re + f_im*f_im);
                fnl_loc[2*i] = fac * f_re;
                fnl_loc[2*i+1] = fac * f_im;
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
              }
            } // ipr

            // array anl_loc is complete

            // compute dfnl[npra][nstloc] = anl^T * c
            if ( basis_.real() )
            {
              // factor 2.0 in next line is: counting G, -G
              // Note: no need to correct for double counting of the
              // G=0 component which is always zero
              double two=1.0;
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
          fion[is][3*ia+0] += tmpfion[3*ia];
          fion[is][3*ia+1] += tmpfion[3*ia+1];
          fion[is][3*ia+2] += tmpfion[3*ia+2];
        }
      }
    } // npr[is]>0
  } // is

  // reduction of enl across rows
  ctxt_.dsum('r',1,1,&enl,1);

  sigma_enl = 0.0;
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
