////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.C,v 1.11 2004-06-01 22:45:30 fgygi Exp $

#include "NonLocalPotential.h"
#include "blas.h"

#if AIX
extern "C" void vsincos(double *x, double *y, double *z, int *n);
#elif OSF1
extern "C" void vcos_sin_(double *x, int *ix, double *y, int *iy, 
              double *z, int *iz, int *n);
#endif


////////////////////////////////////////////////////////////////////////////////
NonLocalPotential::~NonLocalPotential(void)
{
  for ( int is = 0; is < anl.size(); is++ )
  {
    delete anl[is];
  }
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::init(void)
{

  // max deviation from locality at upper quadrature point rcutloc[is]
  const double epsilon = 1.e-4;

  const int ngwl = basis_.localsize();
  
  nsp = atoms_.nsp();
  
  lmax.resize(nsp);
  lloc.resize(nsp);
  lproj.resize(nsp);
  na.resize(nsp);
  naloc.resize(nsp);
  nalocmax.resize(nsp);
  npr.resize(nsp);
  nprna.resize(nsp);
  wt.resize(nsp);
  twnl.resize(nsp);
  dtwnl.resize(nsp);
  
  anl.resize(nsp);
  singr.resize(nsp);
  cosgr.resize(nsp);
  
  nquad.resize(nsp);
  rquad.resize(nsp);
  wquad.resize(nsp);
   
  nspnl = 0;
  for ( int is = 0; is < nsp; is++ )
  {
    Species *s = atoms_.species_list[is];
    
    anl[is] = 0;
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
      
      // block size for distribution of na[is]
      nalocmax[is] = na[is]/ctxt_.npcol() + (na[is]%ctxt_.npcol() > 0 ? 1 : 0);
      
      const int m = 2* sd_.c().m();
      const int mb = 2 * sd_.c().mb();
      const int n = nprna[is];
      const int nb = nalocmax[is] * npr[is];
      anl[is] = new DoubleMatrix(ctxt_,m,n,mb,nb);
      anl[is]->clear(); // clear to zero padding areas
      
      assert(anl[is]->nloc() % npr[is] == 0);
      naloc[is] = anl[is]->nloc() / npr[is];
      singr[is].resize(naloc[is]*ngwl);
      cosgr[is].resize(naloc[is]*ngwl);
      
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
        //cout << " NonLocalPotential::init: trapezoidal rule with right endpoint"
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
    
  const int ngwl = basis_.localsize();
  const double pi = M_PI;
  const double fpi = 4.0 * pi;
  const double s14pi = sqrt(1.0/fpi);
  const double s34pi = sqrt(3.0/fpi);  
  const double s54pi = sqrt(5.0/fpi);
  const double s20pi = sqrt(20.0*pi);
  const double s20pi3 = sqrt(20.0*pi/3.0);
  const double s3 = sqrt(3.0);

  const double *g   = basis_.g_ptr();
  const double *g2  = basis_.g2_ptr();
  const double *gi  = basis_.gi_ptr();
  const double *g_x = basis_.gx_ptr(0);
  const double *g_y = basis_.gx_ptr(1);
  const double *g_z = basis_.gx_ptr(2);
    
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
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              s->dvnlg(0,g[ig],v,dv);

              t0[ig] = s14pi * v;
              
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tmp = gi[ig] * s14pi * dv;
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
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                // I(l=0) = 4 pi j_l(G r) r
                // twnl[is][ipr][l][ig] = 4 pi j_0(Gr_i) r_i Ylm
                // j_0(Gr) * r = sin(Gr) / G
                // Ylm = s14pi
                const double arg = g[ig] * r;
                // Note: for G=0, gi[0] = 0
                
                const double tgx = g_x[ig];
                const double tgy = g_y[ig];
                const double tgz = g_z[ig];
                const double tgi = gi[ig];
                const double tgi2 = tgi * tgi;
                
                const double ts = sin(arg);
                const double tc = cos(arg);
              
                t0[ig] = fpi * s14pi * ts * tgi;
                
                // dtwnl = fpi s14pi G_i G_j / G (r cos(Gr)/G -sin(Gr)/G^2)
                const double tmp = fpi * s14pi * tgi2 * (r*tc - ts*tgi);
                dt0_xx[ig] = tmp * tgx * tgx;
                dt0_yy[ig] = tmp * tgy * tgy;
                dt0_zz[ig] = tmp * tgz * tgz;
                dt0_xy[ig] = tmp * tgx * tgy;
                dt0_yz[ig] = tmp * tgy * tgz;
                dt0_xz[ig] = tmp * tgx * tgz;
              }
            }
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
              const double tg = g[ig];
              s->dvnlg(l,tg,v,dv);
              
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tgi = gi[ig];
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
                const double tg = g[ig];
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
              
                const double tgx = g_x[ig];
                const double tgy = g_y[ig];
                const double tgz = g_z[ig];
                const double tgx2 = tgx * tgx;
                const double tgy2 = tgy * tgy;
                const double tgz2 = tgz * tgz;
 
                const double tgi = gi[ig];
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
              const double tg = g[ig];
              
              s->dvnlg(l,tg,v,dv);
              
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tgi = gi[ig];
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
                
                const double tg = g[ig];
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
              
                const double tgx = g_x[ig];
                const double tgy = g_y[ig];
                const double tgz = g_z[ig];
                const double tgx2 = tgx * tgx;
                const double tgy2 = tgy * tgy;
                const double tgz2 = tgz * tgz;
 
                const double tgi = gi[ig];
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
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::update_eigr(vector<vector<double> >& tau)
{
  // recompute the values of singr[is][ig+ia*ngwl] and cosgr[is][ig+ia*ngwl]
  // using tau[is][j+3*ia] and gx[ig+j*ngwl]
  // must be called whenever atoms move, or the cell size changes

  int ngw = basis_.size();
  int ngwl = basis_.localsize();
  for ( int is = 0; is < nsp; is++ )
  {
    if ( npr[is] > 0 ) // species is non-local
    {
      int nalocis = naloc[is];
      // compute array gr[ig+ia*ngwl]
      
      int len = nalocis * ngwl;
      double *gr = new double[len]; // gr[ig+ia*ngwl]
      int k = 3;
      double mone = -1.0, zero = 0.0;
      char cn='n';
      
      // next line: const cast is ok since dgemm_ does not modify argument
      double* gx = const_cast<double*>(basis_.gx_ptr(0));
      int iafirst = ctxt_.mycol() * nalocmax[is];
      
      dgemm_(&cn,&cn,&ngwl,&nalocis,&k,&mone,
             gx,&ngwl, &tau[is][3*iafirst],&k, &zero,&gr[0],&ngwl);

#if AIX
      vsincos(&singr[is][0],&cosgr[is][0],&gr[0],&len);
#elif OSF1
      int inc1 = 1;
      vcos_sin_(&gr[0],&inc1,&cosgr[is][0],&inc1,&singr[is][0],&inc1,&len);
#else
      for ( int i = 0; i < len; i++ )
      {
        const double arg = gr[i];
        singr[is][i] = sin(arg);
        cosgr[is][i] = cos(arg);
      }
#endif
      delete [] gr;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void NonLocalPotential::update_anl(void)
{
  // update the projector matrices anl[is]
  
  // cmloc_anl: leading dimension of local array of matrix anl
  const int cmloc_anl = sd_.c().mloc();
  const int ngwl = basis_.localsize();

  for ( int is = 0; is < nsp; is++ )   
  {  
    if ( npr[is] > 0 ) // species is is non-local
    {
      for ( int ipr = 0; ipr < npr[is]; ipr++ )
      {
        const int l = lproj[is][ipr];
        
        // twnl[is][ig+ngwl*ipr]
        const double *t = &twnl[is][ngwl*ipr];
        for ( int ia = 0; ia < naloc[is]; ia++ )
        {
          // anl[is][ig+ipra*ngwl]
          // index = ig+cmloc_anl*(ia+nais*ipr), ig=0
          const int ipra = ia+naloc[is]*ipr;
          double *a = anl[is]->valptr(2*(cmloc_anl*ipra));
          const double *c = &cosgr[is][ia*ngwl];
          const double *s = &singr[is][ia*ngwl];
          
          if ( l == 0 )
          {
          
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              const double tt = t[ig];
              // anl[is][ipr][ia][ig].re =
              //   twnl[is][ipr][ig] * cosgr[is][ia][ig]
              *a++ = tt * *c++;
              // anl[is][ipr][ia][ig].im =
              //   twnl[is][ipr][ig] * singr[is][ia][ig]
              *a++ = tt * *s++;
            }
          }
          else if ( l == 1 )
          {
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              const double tt = t[ig];
              /* Next line: -i * eigr */
              /* -i * (a+i*b) = b - i*a */
              *a++ =  tt * *s++;
              *a++ = -tt * *c++;
            }
          }
          else if ( l == 2 )
          {
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              // Next line: (-) sign for -eigr
              const double tt = t[ig];
              *a++ = - tt * *c++;
              *a++ = - tt * *s++;
            }
          }
        } // ia
      } // ipr
    } // if non-local
  } // is
}

////////////////////////////////////////////////////////////////////////////////
double NonLocalPotential::energy(bool compute_hpsi, SlaterDet& dsd, 
    bool compute_forces, vector<vector<double> >& fion, 
    bool compute_stress, valarray<double>& sigma_enl)
{

  const vector<double>& occ = sd_.occ();

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
      // compute fnl
      // block distribution for fnl: same as SlaterDet for nst
      DoubleMatrix fnl(ctxt_,anl[is]->n(),sd_.c().n(),
                       anl[is]->nb(),sd_.c().nb());
      
      const DoubleMatrix c_proxy(sd_.c());
      
      fnl.gemm('t','n',2.0,*anl[is],c_proxy,0.0);
      
      // correct for double counting of G=0 components
      // rank-1 update using first row of *anl[is] and c_proxy
      fnl.ger(-1.0,*anl[is],0,c_proxy,0);
      
      // compute the non-local energy
      // multiply fnl[ipra+nprna*n] by fac = wt[is][ipr] * omega_inv;
      // block sizes: npr*nalocmax x c().nb()
      // loop over local array
      double*f = fnl.valptr(0);
      const int mb = fnl.mb();
      const int nb = fnl.nb();
      const int mloc = fnl.mloc();
      for ( int li=0; li < fnl.mblocks(); li++)
      {
        const int mbs = fnl.mbs(li);
        for ( int lj=0; lj < fnl.nblocks(); lj++)
        {
          const int nbs = fnl.nbs(lj);
          for ( int ii=0; ii < mbs; ii++)
          {
            const int ipr = ii / nalocmax[is];
            const double fac = wt[is][ipr] * omega_inv;
            for ( int jj=0; jj < nbs; jj++)
            {
              // global index: i(li,ii), j(lj,jj)
              const int nglobal = fnl.j(lj,jj);
              const double facn = fac * occ[nglobal];
              const int iii = ii+li*mb;
              const int jjj = jj+lj*nb;
              const double tmp = f[iii+mloc*jjj];
              enl += facn * tmp * tmp;
              f[iii+mloc*jjj] = fac * tmp;
            }
          }
        }
      }

      if ( compute_hpsi )
      {
        // Apply operator to electronic states and accumulate in dsd
        DoubleMatrix cp_proxy(dsd.c());
        cp_proxy.gemm('n','n',1.0,*anl[is],fnl,1.0);
      }

      // ionic forces
      if ( compute_forces )
      {
        double *tmpfion = new double[3*na[is]];
        for ( int i = 0; i < 3*na[is]; i++ )
          tmpfion[i] = 0.0;
          
        DoubleMatrix danl(ctxt_,anl[is]->m(),anl[is]->n(),
          anl[is]->mb(),anl[is]->nb());
        DoubleMatrix dfnl(ctxt_,fnl.m(),fnl.n(),fnl.mb(),fnl.nb());
        const int ngwl = basis_.localsize();
          
        for ( int j = 0; j < 3; j++ )
        {
          const double *const gxj = basis_.gx_ptr(j);
          for ( int ipr = 0; ipr < npr[is]; ipr++ )
          {
            const int l = lproj[is][ipr];
 
            // twnl[is][ig+ngwl*ipr]
            const double *t = &twnl[is][ngwl*ipr];
            for ( int ia = 0; ia < naloc[is]; ia++ )
            {
              // danl[ig+ipra*ngwl]
              // index = ig+cmloc_anl*(ia+nais*ipr), ig=0
              const int ipra = ia+naloc[is]*ipr;
              double *da = danl.valptr(2*(sd_.c().mloc()*ipra));
              const double *c = &cosgr[is][ia*ngwl];
              const double *s = &singr[is][ia*ngwl];
 
              if ( l == 0 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double tt = gxj[ig] * t[ig];
                  // Next lines: -i * ( a + ib ) = b - ia
                  *da++ =  tt * *s++;
                  *da++ = -tt * *c++;
                }
              }
              else if ( l == 1 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  // Next lines: (-i)**2 * ( a + ib ) = - a - ib
                  const double tt = - gxj[ig] * t[ig];
                  *da++ = tt * *c++;
                  *da++ = tt * *s++;
                }
              }
              else if ( l == 2 )
              {
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  // Next lines: (-i) * - ( a + ib ) = i*(a+ib) = - b + ia
                  const double tt = gxj[ig] * t[ig];
                  *da++ = -tt * *s++;
                  *da++ =  tt * *c++;
                }
              }
            } // ia
          } // ipr

          // compute dfnl
          const DoubleMatrix c_proxy(sd_.c());
 
          dfnl.gemm('t','n',2.0,danl,c_proxy,0.0);
 
          // Note: no need to correct for double counting of the
          // G=0 component which is always zero

          // non-local forces
 
          // loop over local array
          // block sizes: npr*nalocmax x c().nb()
          const double*f = fnl.valptr(0);
          const double*df = dfnl.valptr(0);
          const int mloc = fnl.mloc();
          const int mb = fnl.mb();
          const int nb = fnl.nb();
          const int ia_first = fnl.context().myrow() * nalocmax[is];
          for ( int li=0; li < fnl.mblocks(); li++)
          {
            const int mbs = fnl.mbs(li);
            for ( int lj=0; lj < fnl.nblocks(); lj++)
            {
              const int nbs = fnl.nbs(lj);
              for ( int ii=0; ii < mbs; ii++)
              {
                const int ia_local = ii % nalocmax[is];
                // ia_global = ia_local + myrow * nalocmax[is]
                const int ia_global = ia_local + ia_first;
                for ( int jj=0; jj < nbs; jj++)
                {
                  const int nglobal = fnl.j(lj,jj);
                  // Factor 2.0 in next line from derivative of |Fnl|^2
                  const double facn = 2.0 * occ[nglobal];
                  const int iii = ii+li*mb;
                  const int jjj = jj+lj*nb;
                  tmpfion[3*ia_global+j] -= facn *
                    f[iii+mloc*jjj] * df[iii+mloc*jjj];
                }
              }
            }
          }
        } // j
        
        ctxt_.dsum(3*na[is],1,tmpfion,3*na[is]);
        for ( int ia = 0; ia < na[is]; ia++ )
        {
          fion[is][3*ia+0] += tmpfion[3*ia];
          fion[is][3*ia+1] += tmpfion[3*ia+1];
          fion[is][3*ia+2] += tmpfion[3*ia+2];
        }
        delete [] tmpfion;
      } // compute_forces
        
      if ( compute_stress )
      {
        const int ngwl = basis_.localsize();
        DoubleMatrix danl(ctxt_,anl[is]->m(),anl[is]->n(),
          anl[is]->mb(),anl[is]->nb());
        DoubleMatrix dfnl(ctxt_,fnl.m(),fnl.n(),fnl.mb(),fnl.nb());
        
        for ( int ij = 0; ij < 6; ij++ )
        {
          int ipr = 0;
          while ( ipr < npr[is] )
          {
            const int l = lproj[is][ipr];
            if ( l == 0 )
            {
              // dtwnl[is][ipr][ij][ngwl]
              // index = ig + ngwl * ( ij + 6 * ipr))
              // ipr = iquad + nquad[is] * ilm, where ilm = 0
              const double *const dt0 = &dtwnl[is][ngwl*(ij+6*ipr)];
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra0 = ia+naloc[is]*ipr;
                double *da0 = danl.valptr(2*(sd_.c().mloc()*ipra0));
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d0 = dt0[ig];
                  // danl[is][ipr][iquad][ia][ig].re =
                  //   dtwnl[is][ipr][iquad][j][ig] * cosgr[is][ia][ig]
                  *da0++ = *c++ * d0;
                  // danl[is][ipr][iquad][ia][ig].im =
                  //   dtwnl[is][ipr][iquad][j][ig] * singr[is][ia][ig]
                  *da0++ = *s++ * d0;
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
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra1 = ia+naloc[is]*ipr1;
                const int ipra2 = ia+naloc[is]*ipr2;
                const int ipra3 = ia+naloc[is]*ipr3;
                double *da1 = danl.valptr(2*(sd_.c().mloc()*ipra1));
                double *da2 = danl.valptr(2*(sd_.c().mloc()*ipra2));
                double *da3 = danl.valptr(2*(sd_.c().mloc()*ipra3));
 
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d1 = dt1[ig];
                  const double d2 = dt2[ig];
                  const double d3 = dt3[ig];
                  // Next line: (-i)^l factor is -i
                  // Next line: -i * eigr
                  // -i * (a+i*b) = b - i*a
                  const double tc = -*c++; //  -cosgr[is][ia][ig]
                  const double ts =  *s++; //   singr[is][ia][ig]
                  *da1++ = d1 * ts;
                  *da1++ = d1 * tc;
                  *da2++ = d2 * ts;
                  *da2++ = d2 * tc;
                  *da3++ = d3 * ts;
                  *da3++ = d3 * tc;
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
              for ( int ia = 0; ia < naloc[is]; ia++ )
              {
                const int ipra4 = ia+naloc[is]*ipr4;
                const int ipra5 = ia+naloc[is]*ipr5;
                const int ipra6 = ia+naloc[is]*ipr6;
                const int ipra7 = ia+naloc[is]*ipr7;
                const int ipra8 = ia+naloc[is]*ipr8;
                double *da4 = danl.valptr(2*(sd_.c().mloc()*ipra4));
                double *da5 = danl.valptr(2*(sd_.c().mloc()*ipra5));
                double *da6 = danl.valptr(2*(sd_.c().mloc()*ipra6));
                double *da7 = danl.valptr(2*(sd_.c().mloc()*ipra7));
                double *da8 = danl.valptr(2*(sd_.c().mloc()*ipra8));
 
                const double *c = &cosgr[is][ia*ngwl];
                const double *s = &singr[is][ia*ngwl];
                for ( int ig = 0; ig < ngwl; ig++ )
                {
                  const double d4 = dt4[ig];
                  const double d5 = dt5[ig];
                  const double d6 = dt6[ig];
                  const double d7 = dt7[ig];
                  const double d8 = dt8[ig];
                  // Next lines: (-i)^2 * ( a + ib ) =  - ( a + ib )
                  const double tc = -*c++;
                  const double ts = -*s++;
                  *da4++ = d4 * tc;
                  *da4++ = d4 * ts;
                  *da5++ = d5 * tc;
                  *da5++ = d5 * ts;
                  *da6++ = d6 * tc;
                  *da6++ = d6 * ts;
                  *da7++ = d7 * tc;
                  *da7++ = d7 * ts;
                  *da8++ = d8 * tc;
                  *da8++ = d8 * ts;
                }
              }
            }                                                                
            else                                                               
            {                                                                  
              assert(false);                                                   
            } // l
            ipr += 2*l+1;
          } // while ipr
      
          // compute dfnl
          const DoubleMatrix c_proxy(sd_.c());
 
          dfnl.gemm('t','n',2.0,danl,c_proxy,0.0);
 
          // Note: no need to correct for double counting of the
          // G=0 component which is always zero
          
           // partial contributions to the stress sigma_ij
          // Note: fnl was already premultiplied by the factor
          // fac = wt[is][ipr][iquad] * omega_inv;
          const double *const f = fnl.cvalptr(0);
          const double *const df = dfnl.cvalptr(0);
          const int mb = fnl.mb();
          const int nb = fnl.nb();
          const int mloc = fnl.mloc();
          for ( int li=0; li < fnl.mblocks(); li++)
          {
            const int mbs = fnl.mbs(li);
            for ( int lj=0; lj < fnl.nblocks(); lj++)
            {
              const int nbs = fnl.nbs(lj);
              for ( int ii=0; ii < mbs; ii++)
              {
                const int ipr = ii / nalocmax[is];
                for ( int jj=0; jj < nbs; jj++)
                {
                  // global index: i(li,ii), j(lj,jj)
                  const int nglobal = fnl.j(lj,jj);
                  const double facn = 2.0 * occ[nglobal];
                  const int iii = ii+li*mb;
                  const int jjj = jj+lj*nb;
                  const double tmp = f[iii+mloc*jjj];
                  const double dtmp = df[iii+mloc*jjj];
                  tsum[ij] += facn * tmp * dtmp;
                }
              }
            }
          }
        } // ij
      } // compute_stress
    } // npr[is]>0
  } // is

  ctxt_.dsum(1,1,&enl,1);

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

