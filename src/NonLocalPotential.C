////////////////////////////////////////////////////////////////////////////////
//
// NonLocalPotential.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NonLocalPotential.C,v 1.8 2004-01-22 01:19:34 fgygi Exp $

#include "NonLocalPotential.h"
#include "blas.h"
void gauleg(const double x1,const double x2,double *x,double *w,const int n);
void gaujac(double x1, double x2, double *x, double *w, int n, 
  double alf, double bet);
double gammln(double xx);

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
      dtwnl[is].resize(npr[is]*3*ngwl);
      
      // quadrature abcissae and weights
      rquad[is].resize(nquad[is]);
      wquad[is].resize(nquad[is]);
      
      enum quadrature_rule_type { TRAPEZOID, MODIF_TRAPEZOID, 
      TRAPEZOID_WITH_RIGHT_ENDPOINT, 
      SIMPSON, GAUSS_LEGENDRE, GAUSS_JACOBI };
      
      const quadrature_rule_type quad_rule = TRAPEZOID;
      //const quadrature_rule_type quad_rule = MODIF_TRAPEZOID;
      //const quadrature_rule_type quad_rule = TRAPEZOID_WITH_RIGHT_ENDPOINT;
      //const quadrature_rule_type quad_rule = SIMPSON;
      //const quadrature_rule_type quad_rule = GAUSS_LEGENDRE;
      //const quadrature_rule_type quad_rule = GAUSS_JACOBI;

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
      else if ( quad_rule == GAUSS_LEGENDRE )
      {
        // Generate Gauss-Legendre abscissae and weights
        gauleg(0.0,s->rquad(),&rquad[is][0],&wquad[is][0],nquad[is]);
        //cout << " NonLocalPotential::init: Gauss-Legendre quadrature" << endl;
      }
      else if ( quad_rule == GAUSS_JACOBI )
      {
        // Generate Gauss-Jacobi abscissae and weights
        // int (1-x)^alpha (1+x)^beta f(x)
        const double alpha = 0.0;
        const double beta = 1.0;
        gaujac(0.0,s->rquad(),&rquad[is][0],&wquad[is][0],nquad[is],
               alpha, beta);
        //cout << " NonLocalPotential::init: Gauss-Jacobi quadrature" << endl;
        //cout << " alpha = " << alpha << " beta = " << beta << endl;
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
            // dtwnl[is][ipr][j][ngwl]
            // index = ig + ngwl * ( j + 3 * ipr ), ipr = 0
            // i.e. index = ig + j * ngwl
            double *dt0_x = &dtwnl[is][0*ngwl];
            double *dt0_y = &dtwnl[is][1*ngwl];
            double *dt0_z = &dtwnl[is][2*ngwl];
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
              
              const double tmp = - gi[ig] * s14pi * dv;
              dt0_x[ig] = tgx2 * tmp;
              dt0_y[ig] = tgy2 * tmp;
              dt0_z[ig] = tgz2 * tmp;
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
              // index = ig + ngwl * ( j + 3 * iquad)
              double *dt0_x = &dtwnl[is][ngwl*(0+3*iquad)];
              double *dt0_y = &dtwnl[is][ngwl*(1+3*iquad)];
              double *dt0_z = &dtwnl[is][ngwl*(2+3*iquad)];
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
                const double tgx2 = tgx * tgx;
                const double tgy2 = tgy * tgy;
                const double tgz2 = tgz * tgz;
                const double tgi = gi[ig];
                const double tgi2 = tgi * tgi;
                
                const double ts = sin(arg);
                const double tc = cos(arg);
              
                t0[ig] = fpi * s14pi * ts * tgi;
                
                // dtwnl = - fpi s14pi G_j^2 / G (r cos(Gr)/G -sin(Gr)/G^2)
                const double tmp = - fpi * s14pi * tgi2 * (r*tc - ts*tgi);
                dt0_x[ig] = tgx2 * tmp;
                dt0_y[ig] = tgy2 * tmp;
                dt0_z[ig] = tgz2 * tmp;
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
            
            // dtwnl[is][ipr][j][ngwl]
            // index = ig + ngwl * ( j + 3 * ipr )
            double *dt1_x = &dtwnl[is][ngwl*(0+3*ipr1)];
            double *dt1_y = &dtwnl[is][ngwl*(1+3*ipr1)];
            double *dt1_z = &dtwnl[is][ngwl*(2+3*ipr1)];

            double *dt2_x = &dtwnl[is][ngwl*(0+3*ipr2)];
            double *dt2_y = &dtwnl[is][ngwl*(1+3*ipr2)];
            double *dt2_z = &dtwnl[is][ngwl*(2+3*ipr2)];

            double *dt3_x = &dtwnl[is][ngwl*(0+3*ipr3)];
            double *dt3_y = &dtwnl[is][ngwl*(1+3*ipr3)];
            double *dt3_z = &dtwnl[is][ngwl*(2+3*ipr3)];
            
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              s->dvnlg(l,g[ig],v,dv);
              
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

              // dylm/dx, dylm/dy, c * dylm/dz
              const double dy1_x = s34pi * tgi * ( 1.0 - tgx2 * tgi2 );
              const double dy1_y = s34pi * tgi * ( - tgx * tgy * tgi2 );
              const double dy1_z = s34pi * tgi * ( - tgx * tgz * tgi2 );
              
              const double dy2_x = s34pi * tgi * ( - tgy * tgx * tgi2 );
              const double dy2_y = s34pi * tgi * ( 1.0 - tgy2 * tgi2 );
              const double dy2_z = s34pi * tgi * ( - tgy * tgz * tgi2 );
              
              const double dy3_x = s34pi * tgi * ( - tgz * tgx * tgi2 );
              const double dy3_y = s34pi * tgi * ( - tgz * tgy * tgi2 );
              const double dy3_z = s34pi * tgi * ( 1.0 - tgz2 * tgi2 );
              
              dt1_x[ig] = - ( v * tgx * dy1_x + y1 * tgx2 * tgi * dv );
              dt1_y[ig] = - ( v * tgy * dy1_y + y1 * tgy2 * tgi * dv );
              dt1_z[ig] = - ( v * tgz * dy1_z + y1 * tgz2 * tgi * dv );

              dt2_x[ig] = - ( v * tgx * dy2_x + y2 * tgx2 * tgi * dv );
              dt2_y[ig] = - ( v * tgy * dy2_y + y2 * tgy2 * tgi * dv );
              dt2_z[ig] = - ( v * tgz * dy2_z + y2 * tgz2 * tgi * dv );

              dt3_x[ig] = - ( v * tgx * dy3_x + y3 * tgx2 * tgi * dv );
              dt3_y[ig] = - ( v * tgy * dy3_y + y3 * tgy2 * tgi * dv );
              dt3_z[ig] = - ( v * tgz * dy3_z + y3 * tgz2 * tgi * dv );
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
              // index = ig + ngwl * ( j + 3 * ( iquad + nquad[is] * ilm ))
              double *dt1_x = &dtwnl[is][ngwl*(0+3*ipr1)];
              double *dt1_y = &dtwnl[is][ngwl*(1+3*ipr1)];
              double *dt1_z = &dtwnl[is][ngwl*(2+3*ipr1)];

              double *dt2_x = &dtwnl[is][ngwl*(0+3*ipr2)];
              double *dt2_y = &dtwnl[is][ngwl*(1+3*ipr2)];
              double *dt2_z = &dtwnl[is][ngwl*(2+3*ipr2)];

              double *dt3_x = &dtwnl[is][ngwl*(0+3*ipr3)];
              double *dt3_y = &dtwnl[is][ngwl*(1+3*ipr3)];
              double *dt3_z = &dtwnl[is][ngwl*(2+3*ipr3)];
            
              const double r = rquad[is][iquad];
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                double v = 0.0, dv = 0.0;
                // j_1(Gr) = (1/(Gr))*(sin(Gr)/(Gr)-cos(Gr))
                const double z = g[ig] * r;
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

                // dylm/dx, dylm/dy, c * dylm/dz
                const double dy1_x = s34pi * tgi * ( 1.0 - tgx2 * tgi2 );
                const double dy1_y = s34pi * tgi * ( - tgx * tgy * tgi2 );
                const double dy1_z = s34pi * tgi * ( - tgx * tgz * tgi2 );
 
                const double dy2_x = s34pi * tgi * ( - tgy * tgx * tgi2 );
                const double dy2_y = s34pi * tgi * ( 1.0 - tgy2 * tgi2 );
                const double dy2_z = s34pi * tgi * ( - tgy * tgz * tgi2 );
 
                const double dy3_x = s34pi * tgi * ( - tgz * tgx * tgi2 );
                const double dy3_y = s34pi * tgi * ( - tgz * tgy * tgi2 );
                const double dy3_z = s34pi * tgi * ( 1.0 - tgz2 * tgi2 );
 
                dt1_x[ig] = - ( v * tgx * dy1_x + y1 * tgx2 * tgi * dv );
                dt1_y[ig] = - ( v * tgy * dy1_y + y1 * tgy2 * tgi * dv );
                dt1_z[ig] = - ( v * tgz * dy1_z + y1 * tgz2 * tgi * dv );

                dt2_x[ig] = - ( v * tgx * dy2_x + y2 * tgx2 * tgi * dv );
                dt2_y[ig] = - ( v * tgy * dy2_y + y2 * tgy2 * tgi * dv );
                dt2_z[ig] = - ( v * tgz * dy2_z + y2 * tgz2 * tgi * dv );

                dt3_x[ig] = - ( v * tgx * dy3_x + y3 * tgx2 * tgi * dv );
                dt3_y[ig] = - ( v * tgy * dy3_y + y3 * tgy2 * tgi * dv );
                dt3_z[ig] = - ( v * tgz * dy3_z + y3 * tgz2 * tgi * dv );
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
            
            // dtwnl[is][ipr][j][ngwl]
            // index = ig + ngwl * ( j + 3 * ipr )
            double *dt4_x = &dtwnl[is][ngwl*(0+3*ipr4)];
            double *dt4_y = &dtwnl[is][ngwl*(1+3*ipr4)];
            double *dt4_z = &dtwnl[is][ngwl*(2+3*ipr4)];

            double *dt5_x = &dtwnl[is][ngwl*(0+3*ipr5)];
            double *dt5_y = &dtwnl[is][ngwl*(1+3*ipr5)];
            double *dt5_z = &dtwnl[is][ngwl*(2+3*ipr5)];

            double *dt6_x = &dtwnl[is][ngwl*(0+3*ipr6)];
            double *dt6_y = &dtwnl[is][ngwl*(1+3*ipr6)];
            double *dt6_z = &dtwnl[is][ngwl*(2+3*ipr6)];
            
            double *dt7_x = &dtwnl[is][ngwl*(0+3*ipr7)];
            double *dt7_y = &dtwnl[is][ngwl*(1+3*ipr7)];
            double *dt7_z = &dtwnl[is][ngwl*(2+3*ipr7)];
            
            double *dt8_x = &dtwnl[is][ngwl*(0+3*ipr8)];
            double *dt8_y = &dtwnl[is][ngwl*(1+3*ipr8)];
            double *dt8_z = &dtwnl[is][ngwl*(2+3*ipr8)];
            
            for ( int ig = 0; ig < ngwl; ig++ )
            {
              double v,dv;
              s->dvnlg(l,g[ig],v,dv);
              
              const double tgx = g_x[ig];
              const double tgy = g_y[ig];
              const double tgz = g_z[ig];
              const double tgx2 = tgx * tgx;
              const double tgy2 = tgy * tgy;
              const double tgz2 = tgz * tgz;
              
              const double tgi = gi[ig];
              const double tg2 = g2[ig];
              const double tgi2 = tgi * tgi;
              
              const double s3 = sqrt(3.0);
              const double yfac = s54pi * tgi2;
              
              const double y4 = 0.5 * (3.0 * tgz2 - tg2 ) * yfac;
              const double y5 = 0.5 * s3 * ( tgx2 - tgy2 ) * yfac;
              const double y6 = s3 * tgx * tgy * yfac;
              const double y7 = s3 * tgy * tgz * yfac;
              const double y8 = s3 * tgz * tgx * yfac;
              
              // dylm/dx, dylm/dy, c * dylm/dz
              
              const double tgi4 = tgi2 * tgi2;
              
              // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
              const double dy4_x = - s54pi *3.0 * tgx * tgz2 * tgi4;
              const double dy4_y = - s54pi *3.0 * tgy * tgz2 * tgi4;
              const double dy4_z = - s54pi *3.0 * tgz * ( tgz2 - tg2 ) * tgi4;
              
              // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
              const double dy5_x = - s54pi *s3 * tgx * ( tgx2 - tgy2 - tg2 ) * tgi4;
              const double dy5_y = - s54pi *s3 * tgy * ( tgx2 - tgy2 + tg2 ) * tgi4;
              const double dy5_z = - s54pi *s3 * tgz * ( tgx2 - tgy2 ) * tgi4;
              
              // y6 = s54pi sqrt(3) x y / r^2
              const double dy6_x = - s54pi *s3 * tgy * ( 2 * tgx2 - tg2 ) * tgi4;
              const double dy6_y = - s54pi *s3 * tgx * ( 2 * tgy2 - tg2 ) * tgi4;
              const double dy6_z = - s54pi *s3 * tgz * 2 * tgx * tgy * tgi4;
              
              // y7 = s54pi sqrt(3) y z / r^2
              const double dy7_x = - s54pi *s3 * tgx * 2 * tgy * tgz * tgi4;
              const double dy7_y = - s54pi *s3 * tgz * ( 2 * tgy2 - tg2 ) * tgi4;
              const double dy7_z = - s54pi *s3 * tgy * ( 2 * tgz2 - tg2 ) * tgi4;
              
              // y8 = s54pi sqrt(3) z x / r^2
              const double dy8_x = - s54pi *s3 * tgz * ( 2 * tgx2 - tg2 ) * tgi4;
              const double dy8_y = - s54pi *s3 * tgy * 2 * tgz * tgx * tgi4;
              const double dy8_z = - s54pi *s3 * tgx * ( 2 * tgz2 - tg2 ) * tgi4;
              
              t4[ig]  = y4 * v;
              t5[ig]  = y5 * v;
              t6[ig]  = y6 * v;
              t7[ig]  = y7 * v;
              t8[ig]  = y8 * v;

              dt4_x[ig] = - ( v * tgx * dy4_x + y4 * tgx2 * tgi * dv );
              dt4_y[ig] = - ( v * tgy * dy4_y + y4 * tgy2 * tgi * dv );
              dt4_z[ig] = - ( v * tgz * dy4_z + y4 * tgz2 * tgi * dv );

              dt5_x[ig] = - ( v * tgx * dy5_x + y5 * tgx2 * tgi * dv );
              dt5_y[ig] = - ( v * tgy * dy5_y + y5 * tgy2 * tgi * dv );
              dt5_z[ig] = - ( v * tgz * dy5_z + y5 * tgz2 * tgi * dv );

              dt6_x[ig] = - ( v * tgx * dy6_x + y6 * tgx2 * tgi * dv );
              dt6_y[ig] = - ( v * tgy * dy6_y + y6 * tgy2 * tgi * dv );
              dt6_z[ig] = - ( v * tgz * dy6_z + y6 * tgz2 * tgi * dv );

              dt7_x[ig] = - ( v * tgx * dy7_x + y7 * tgx2 * tgi * dv );
              dt7_y[ig] = - ( v * tgy * dy7_y + y7 * tgy2 * tgi * dv );
              dt7_z[ig] = - ( v * tgz * dy7_z + y7 * tgz2 * tgi * dv );

              dt8_x[ig] = - ( v * tgx * dy8_x + y8 * tgx2 * tgi * dv );
              dt8_y[ig] = - ( v * tgy * dy8_y + y8 * tgy2 * tgi * dv );
              dt8_z[ig] = - ( v * tgz * dy8_z + y8 * tgz2 * tgi * dv );
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

              // dtwnl[is][ipr][j][ngwl]
              // index = ig + ngwl * ( j + 3 * ( iquad + nquad[is] * ipr ))
              double *dt4_x = &dtwnl[is][ngwl*(0+3*ipr4)];
              double *dt4_y = &dtwnl[is][ngwl*(1+3*ipr4)];
              double *dt4_z = &dtwnl[is][ngwl*(2+3*ipr4)];

              double *dt5_x = &dtwnl[is][ngwl*(0+3*ipr5)];
              double *dt5_y = &dtwnl[is][ngwl*(1+3*ipr5)];
              double *dt5_z = &dtwnl[is][ngwl*(2+3*ipr5)];

              double *dt6_x = &dtwnl[is][ngwl*(0+3*ipr6)];
              double *dt6_y = &dtwnl[is][ngwl*(1+3*ipr6)];
              double *dt6_z = &dtwnl[is][ngwl*(2+3*ipr6)];

              double *dt7_x = &dtwnl[is][ngwl*(0+3*ipr7)];
              double *dt7_y = &dtwnl[is][ngwl*(1+3*ipr7)];
              double *dt7_z = &dtwnl[is][ngwl*(2+3*ipr7)];

              double *dt8_x = &dtwnl[is][ngwl*(0+3*ipr8)];
              double *dt8_y = &dtwnl[is][ngwl*(1+3*ipr8)];
              double *dt8_z = &dtwnl[is][ngwl*(2+3*ipr8)];
            
              const double r = rquad[is][iquad];
              for ( int ig = 0; ig < ngwl; ig++ )
              {
                double v = 0.0, dv = 0.0;
                // j_2(z) = (3/z^3-1/z) sin(z) - 3/z^2 cos(z)
                // j_2(z) = (1/z)*((3/z^2-1)*sin(z) - (3/z) cos(z) )
                const double z = g[ig] * r;
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
                const double tg2 = g2[ig];
                const double tgi2 = tgi * tgi;
 
                const double s3 = sqrt(3.0);
                const double yfac = s54pi * tgi2;
 
                const double y4 = 0.5 * (3.0 * tgz2 - tg2 ) * yfac;
                const double y5 = 0.5 * s3 * ( tgx2 - tgy2 ) * yfac;
                const double y6 = s3 * tgx * tgy * yfac;
                const double y7 = s3 * tgy * tgz * yfac;
                const double y8 = s3 * tgz * tgx * yfac;
 
                t4[ig]  = y4 * v;
                t5[ig]  = y5 * v;
                t6[ig]  = y6 * v;
                t7[ig]  = y7 * v;
                t8[ig]  = y8 * v;

                // dylm/dx, dylm/dy, c * dylm/dz
 
                const double tgi4 = tgi2 * tgi2;
 
                // y4 = s54pi 1/2 ( 3 z^2/r^2 - 1 )
                const double dy4_x = - s54pi *3.0 * tgx * tgz2 * tgi4;
                const double dy4_y = - s54pi *3.0 * tgy * tgz2 * tgi4;
                const double dy4_z = - s54pi *3.0 * tgz * ( tgz2 - tg2 ) * tgi4;
 
                // y5 = s54pi sqrt(3)/2 ( x^2 - y^2 ) / r^2
                const double dy5_x = - s54pi *s3 * tgx * ( tgx2 - tgy2 - tg2 ) * tgi4;
                const double dy5_y = - s54pi *s3 * tgy * ( tgx2 - tgy2 + tg2 ) * tgi4;
                const double dy5_z = - s54pi *s3 * tgz * ( tgx2 - tgy2 ) * tgi4;
 
                // y6 = s54pi sqrt(3) x y / r^2
                const double dy6_x = - s54pi *s3 * tgy * ( 2 * tgx2 - tg2 ) * tgi4;
                const double dy6_y = - s54pi *s3 * tgx * ( 2 * tgy2 - tg2 ) * tgi4;
                const double dy6_z = - s54pi *s3 * tgz * 2 * tgx * tgy * tgi4;
 
                // y7 = s54pi sqrt(3) y z / r^2
                const double dy7_x = - s54pi *s3 * tgx * 2 * tgy * tgz * tgi4;
                const double dy7_y = - s54pi *s3 * tgz * ( 2 * tgy2 - tg2 ) * tgi4;
                const double dy7_z = - s54pi *s3 * tgy * ( 2 * tgz2 - tg2 ) * tgi4;
 
                // y8 = s54pi sqrt(3) z x / r^2
                const double dy8_x = - s54pi *s3 * tgz * ( 2 * tgx2 - tg2 ) * tgi4;
                const double dy8_y = - s54pi *s3 * tgy * 2 * tgz * tgx * tgi4;
                const double dy8_z = - s54pi *s3 * tgx * ( 2 * tgz2 - tg2 ) * tgi4;
 
                dt4_x[ig] = - ( v * tgx * dy4_x + y4 * tgx2 * tgi * dv );
                dt4_y[ig] = - ( v * tgy * dy4_y + y4 * tgy2 * tgi * dv );
                dt4_z[ig] = - ( v * tgz * dy4_z + y4 * tgz2 * tgi * dv );

                dt5_x[ig] = - ( v * tgx * dy5_x + y5 * tgx2 * tgi * dv );
                dt5_y[ig] = - ( v * tgy * dy5_y + y5 * tgy2 * tgi * dv );
                dt5_z[ig] = - ( v * tgz * dy5_z + y5 * tgz2 * tgi * dv );

                dt6_x[ig] = - ( v * tgx * dy6_x + y6 * tgx2 * tgi * dv );
                dt6_y[ig] = - ( v * tgy * dy6_y + y6 * tgy2 * tgi * dv );
                dt6_z[ig] = - ( v * tgz * dy6_z + y6 * tgz2 * tgi * dv );

                dt7_x[ig] = - ( v * tgx * dy7_x + y7 * tgx2 * tgi * dv );
                dt7_y[ig] = - ( v * tgy * dy7_y + y7 * tgy2 * tgi * dv );
                dt7_z[ig] = - ( v * tgz * dy7_z + y7 * tgz2 * tgi * dv );

                dt8_x[ig] = - ( v * tgx * dy8_x + y8 * tgx2 * tgi * dv );
                dt8_y[ig] = - ( v * tgy * dy8_y + y8 * tgy2 * tgi * dv );
                dt8_z[ig] = - ( v * tgz * dy8_z + y8 * tgz2 * tgi * dv );
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
    bool compute_stress, UnitCell& dcell)
{

  assert(compute_stress==false); // stress for arbitrary cells not implemented
  
  const vector<double>& occ = sd_.occ();

  double enl = 0.0;

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
          
        for ( int j = 0; j < 3; j++ )
        {
          const double *gxj = basis_.gx_ptr(j);
 
          const int ngwl = basis_.localsize();

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
        
#if 0
      if ( ttstress )
      {
        for ( int j = 0; j < 3; j++ )
        { 
          int ipr = 0;
          for ( int l = 0; l <= lmax[is]; l++ )
          {
            if ( l != lloc[is] )
            {
              if ( l == 0 )
              {
                for ( int iquad = 0; iquad < nquad[is]; iquad++ )
                {
                  // dtwnl[is][ipr][iquad][j][ngwl]
                  // index = ig + ngwl * ( j + 3 * ( iquad + nquad[is] * ipr ))
                  const double *dt0 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr))];
                  for ( int ia = 0; ia < nais; ia++ )
                  {
                    // anl[is][ipr][iquad][ia][ig]
                    // index = 
                    // ig + ngwl * ( ia + na[is] * ( iquad + nquad * ipr ))
                    double *a = 
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr)))];
                    const double *c = &cosgr[is][ia*ngwl];
                    const double *s = &singr[is][ia*ngwl];
                    for ( int ig = 0; ig < ngwl; ig++ )
                    {
                      const double d0 = dt0[ig];
                      // anl[is][ipr][iquad][ia][ig].re =
                      //   dtwnl[is][ipr][iquad][j][ig] * cosgr[is][ia][ig]
                      *a++ = *c++ * d0;
                      // anl[is][ipr][iquad][ia][ig].im =
                      //   dtwnl[is][ipr][iquad][j][ig] * singr[is][ia][ig]
                      *a++ = *s++ * d0;
                    }
                  }
                }
                ipr += (2*l+1);
              }
              else if ( l == 1 )
              {
                const int ipr1 = ipr;
                const int ipr2 = ipr + 1;
                const int ipr3 = ipr + 2;
                for ( int iquad = 0; iquad < nquad[is]; iquad++ )
                {
                  // dtwnl[is][ipr][iquad][j][ngwl]
                  // index = ig + ngwl * ( j + 3 * ( iquad + nquad[is] * ipr ))
                  const double *dt1 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr1))];
                  const double *dt2 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr2))];
                  const double *dt3 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr3))];
                  for ( int ia = 0; ia < nais; ia++ )
                  {
                    double *a1 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr1)))];
                    double *a2 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr2)))];
                    double *a3 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr3)))];
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
                      *a1++ = d1 * ts;
                      *a1++ = d1 * tc;
                      *a2++ = d2 * ts;
                      *a2++ = d2 * tc;
                      *a3++ = d3 * ts;
                      *a3++ = d3 * tc;
                    }
                  }
                }
                ipr += (2*l+1);
              }
              else if ( l == 2 )
              {
                const int ipr4 = ipr;
                const int ipr5 = ipr + 1;
                const int ipr6 = ipr + 2;
                const int ipr7 = ipr + 3;
                const int ipr8 = ipr + 4;
                for ( int iquad = 0; iquad < nquad[is]; iquad++ )
                {
                  // dtwnl[is][ipr][iquad][j][ngwl]
                  // index = ig + ngwl * ( j + 3 * ( iquad + nquad[is] * ipr ))
                  const double *dt4 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr4))];
                  const double *dt5 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr5))];
                  const double *dt6 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr6))];
                  const double *dt7 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr7))];
                  const double *dt8 = 
                    &dtwnl[is][ngwl*(j+3*(iquad+nquad[is]*ipr8))];
                  for ( int ia = 0; ia < nais; ia++ )
                  {
                    // anl[is][ipr][iquad][ia][ig]
                    double *a4 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr4)))];
                    double *a5 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr5)))];
                    double *a6 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr6)))];
                    double *a7 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr7)))];
                    double *a8 =
                      &anl[is][2*(ngwl*(ia+na[is]*(iquad+nquad[is]*ipr8)))];
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
                      *a4++ = d4 * tc;
                      *a4++ = d4 * ts;
                      *a5++ = d5 * tc;
                      *a5++ = d5 * ts;
                      *a6++ = d6 * tc;
                      *a6++ = d6 * ts;
                      *a7++ = d7 * tc;
                      *a7++ = d7 * ts;
                      *a8++ = d8 * tc;
                      *a8++ = d8 * ts;
                    }
                  }
                }
              }
              else
              {
                assert(false);
              }
            } // l != lloc[is]
          } // l
      
          for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
          {
            int nst = wf.f0[ispin].nst();
            int nproj = npr[is]*nquad[is]*na[is];
            dgemm_(&ct,&cn,&nst,&nproj,&twongwl,&two,
                  (double *) &wf.f0[ispin].c[0],&twongwl,
                  &anl[is][0],&twongwl,
                  &zero,&dfnl[ispin][0],&nst);
                  
            // Note: no need to correct for double counting of the  
            // G=0 component which is always zero
 
            // partial contributions to the ionic forces
            // Note: fnl was already premultiplied by the factor 
            // fac = wt[is][ipr][iquad] * omega_inv;
            int i = 0;
            for ( int ipr = 0; ipr < npr[is]; ipr++ )
            {
              for ( int iquad = 0; iquad < nquad[is]; iquad++ )
              {
                for ( int ia = 0; ia < nais; ia++ )
                {
                  for ( int n = 0; n < nst; n++ )
                  {
                    double facn = 2.0 * wf.occ[ispin][n];
                    // fnl[ispin][ipr][iquad][ia][n]
                    // i = n+nst*(ia+nais*(iquad+nquad*(ipr)))
                    denlda[j] += facn * fnl[ispin][i] * dfnl[ispin][i];
                    i++;
                  }
                }
              }
            }
          } // ispin
        } // j
      } // compute_stress
#endif
    } // npr[is]>0
  } // is

  ctxt_.dsum(1,1,&enl,1);

#if 0
  if ( compute_stress )
  {
#if USE_MPI
    double tsum[3];
    int status = MPI_Allreduce(denlda,&tsum,3,
      MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    assert(status==0);
    denlda[0] = tsum[0];
    denlda[1] = tsum[1];
    denlda[2] = tsum[2];
#endif

    denlda[0] = ( - enl + denlda[0] ) / basis_.cell().x;
    denlda[1] = ( - enl + denlda[1] ) / basis_.cell().y;
    denlda[2] = ( - enl + denlda[2] ) / basis_.cell().z;
  }
#endif

  return enl;
}

////////////////////////////////////////////////////////////////////////////////
void gauleg(const double x1, const double x2, 
  double *x, double *w, const int n)
{
  const double pi = 3.1415926535897932384;
  const double eps = 1.e-14;

  int m = ( n + 1 ) / 2;
  
  const double xm = 0.5 * ( x2 + x1 );
  const double xl = 0.5 * ( x2 - x1 );
  
  for ( int i = 0; i < m; i++ )
  {
    double z = cos ( pi * ( i + 0.75 ) / ( n + 0.5 ) );
    double pp, z1;
    
    do
    {
      double p1 = 1.0;
      double p2 = 0.0;
      double p3;

      for ( int j = 0; j < n; j++ )
      {
        p3 = p2;
        p2 = p1;
        p1 = ( ( 2.0 * (j+1) - 1.0 ) * z * p2 - j * p3 ) / ( j + 1 );
      }
      
      pp = n * ( z * p1 - p2 ) / ( z * z - 1.0 );
      z1 = z;
      z = z1 - p1 / pp;
    } while ( fabs(z-z1) > eps );
    
    x[i] = xm - xl * z;
    x[n-i-1] = xm + xl * z;
    w[i] = 2.0 * xl / ( ( 1.0 - z * z ) * pp * pp );
    w[n-i-1] = w[i];
    
  }
}

////////////////////////////////////////////////////////////////////////////////
void gaujac(double x1, double x2, double *x, double *w, int n, 
  double alf, double bet)
{
  const double eps = 1.e-14;
  const int maxit = 20;
  double gammln(double xx);
  double pp,z,z1;
  double xm = 0.5 * ( x2 + x1 );
  double xl = 0.5 * ( x2 - x1 );
  for ( int i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      double an = alf / n;
      double bn = bet / n;
      double r1 = (1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
      double r2 = 1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
      z = 1.0 - r1/r2;
    }
    else if ( i == 1 )
    {
      double r1 = (4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
      double r2 = 1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
      double r3 = 1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
      z -= (1.0-z)*r1*r2*r3;
    }
    else if ( i == 2 )
    {
      double r1 = (1.67+0.28*alf)/(1.0+0.37*alf);
      double r2 = 1.0+0.22*(n-8.0)/n;
      double r3 = 1.0+8.0*bet/((6.28+bet)*n*n);
      z -= (x[0]-z)*r1*r2*r3;
    }
    else if ( i == n-2 )
    {
      double r1 = (1.0+0.235*bet)/(0.766+0.119*bet);
      double r2 = 1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
      double r3 = 1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
      z += (z-x[n-4])*r1*r2*r3;
    }
    else if ( i == n-1 )
    {
      double r1 = (1.0+0.37*bet)/(1.67+0.28*bet);
      double r2 = 1.0/(1.0+0.22*(n-8.0)/n);
      double r3 = 1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
      z += (z-x[n-3])*r1*r2*r3;
    }
    else
    {
      z = 3.0 * x[i-1] - 3.0 * x[i-2] + x[i-3];
    }
    
    double alfbet = alf + bet;
    bool done = false;
    double temp, p1, p2, p3;
    for ( int its = 0; its < maxit && !done; its++ )
    {
      temp = 2.0 + alfbet;
      p1 = ( alf - bet + temp * z ) / 2.0;
      p2 = 1.0;
      p3;
      for ( int j = 2; j <= n; j++ )
      {
        p3 = p2;
        p2 = p1;
        temp = 2 * j + alfbet;
        double a = 2 * j * ( j + alfbet ) * ( temp - 2.0 );
        double b = ( temp - 1.0 ) * ( alf*alf - bet*bet + temp * 
                   ( temp - 2.0 ) * z );
        double c = 2.0 * ( j - 1 + alf ) * ( j - 1 + bet ) * temp;
        p1 = ( b * p2 - c * p3 ) / a;
      }
      
      // p1 is now the Jacobi polynomial
      pp = ( n * ( alf - bet - temp * z ) * p1 +
           2.0 * ( n + alf ) * ( n + bet ) * p2 ) / ( temp * ( 1.0 - z*z ) );
      z1 = z;
      z = z1 - p1 / pp;
      
      if ( fabs(z-z1) <= eps )
        done = true;
    }
    
    if ( !done )
    {
      cout << " gauss-jacobi: error, more than " << maxit << " iterations"
      << endl;
      exit(1);
    }
    
    x[i] = z;
    w[i] = exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
               gammln(n+alfbet+1.0)) * temp*pow(2.0,alfbet) / ( pp*p2 );
               
  }
  for ( int i = 0; i < n; i++ )
  {
    // translate nodes to the [x1,x2]
    x[i] = xm - xl * x[i];
    w[i] = xl * w[i];
  }
}

////////////////////////////////////////////////////////////////////////////////
double gammln(double xx)
{
  static double cof[6] = 
  {
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  };
  double y = xx;
  double x = xx;
  double tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  double ser = 1.000000000190015;
  for ( int j = 0; j <= 5; j++ ) ser += cof[j]/++y;
  
  return -tmp+log(2.5066282746310005*ser/x);
}
