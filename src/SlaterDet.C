////////////////////////////////////////////////////////////////////////////////
//
// SlaterDet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SlaterDet.C,v 1.20 2003-11-21 19:12:09 fgygi Exp $

#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Context.h"
#include "blas.h" // daxpy

#include <cstdlib>
#include <iostream>
#include <iomanip>
#if USE_CSTDIO_LFS
#include <sstream>
#include <cstdio>
#endif
using namespace std;

// XML transcoding for print function
#include <xercesc/util/Base64.hpp>
#include <xercesc/util/XMLString.hpp>
using namespace xercesc;

////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(const Context& ctxt, D3vector kpoint) : ctxt_(ctxt), 
basis_(Context(ctxt,'c',ctxt.mycol()),kpoint), c_(ctxt)
{}

////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(const SlaterDet& rhs) : ctxt_(rhs.context()),
basis_(rhs.basis_), c_(rhs.c_) {}
  
////////////////////////////////////////////////////////////////////////////////
SlaterDet::~SlaterDet() 
{ 
  // cout << ctxt_.mype() << ": SlaterDet::dtor: ctxt=" << ctxt_;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::byteswap_double(size_t n, double* x)
{
  if (n==0) return;
  unsigned char* c = (unsigned char*) x;
  while ( n-- > 0 )
  {
    unsigned char tmp;
    tmp = c[7]; c[7] = c[0]; c[0] = tmp;
    tmp = c[6]; c[6] = c[1]; c[1] = tmp;
    tmp = c[5]; c[5] = c[2]; c[2] = tmp;
    tmp = c[4]; c[4] = c[3]; c[3] = tmp;
    
    c+=8;
  }
}  

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::resize(const UnitCell& cell, const UnitCell& refcell, 
  double ecut, int nst)
{
  try
  {
    tmap["resize"].start();
    basis_.resize(cell,refcell,ecut);

    occ_.resize(nst);
    eig_.resize(nst);
    
    int cm = ctxt_.nprow() * basis_.maxlocalsize();
    int nb = nst/ctxt_.npcol() + (nst%ctxt_.npcol() > 0 ? 1 : 0);
    c_.resize(cm,nst,basis_.maxlocalsize(),nb);
    
    if ( nst <= basis_.size() )
    {
      // initialize c_
      c_.clear();
      const double s2i = 1.0 / sqrt(2.0);
      
      // for each n, find the smallest g vector and initialize
      int ismallest = 0;
      // on each process, basis.isort(ismallest) is the index of the smallest
      // local g vector
      for ( int n = 0; n < nst; n++ )
      {
        double value = 1.0;
        if ( basis().real() && n != 0 )
          value = s2i;
          
        // find process row holding the smallest g vector
        double g2 = basis_.g2(basis_.isort(ismallest));
        // cout << "smallest vector on proc " << ctxt_.mype()
        //      << " has norm " << g2 << endl;
        int minrow, mincol;
        ctxt_.dmin('c',' ',1,1,&g2,1,&minrow,&mincol,1,-1,-1);
        
        // find column hosting state n
        int pc = c_.pc(n);
        int pr = minrow;
        if ( pr == ctxt_.myrow() )
        {
          int iii = basis_.isort(ismallest);
          ismallest++; // increment on entire process row
          if ( pc == ctxt_.mycol() )
          {
            // cout << " n=" << n << " on process " 
            //      << pr << "," << pc
            //      << " vector " << basis_.idx(3*iii) << " "
            //      << basis_.idx(3*iii+1) << " "
            //      << basis_.idx(3*iii+2) << " norm=" 
            //      << basis_.g2(iii) << " "
            //      << value << endl;
            int jjj = c_.m(n) * c_.nb() + c_.y(n);
            int index = iii+c_.mloc()*jjj;
            c_[index] = complex<double> (value,0.0);
          }
        }
        
      }
        
    }
    tmap["resize"].stop();
  }
  catch ( bad_alloc )
  {
    cout << " bad_alloc exception caught in SlaterDet::resize" << endl;
    throw;
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::ekin(void) const
{
  double sum = 0.0;
  const double* kpg2 = basis_.kpg2_ptr();
  
  // factor in next line: 2.0 for G and -G (if basis is real) and
  // 0.5 from 1/(2m)
  // note: if basis is real, the factor of 2.0 for G=0 need not be
  // corrected since G^2 = 0
  
  const double fac = basis_.real() ? 1.0 : 0.5;
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    // global n index
    int nn = ctxt_.mycol() * c_.nb() + n;
    
    if ( occ_[nn] != 0.0 )
    {
      const complex<double>* p = c_.cvalptr(c_.mloc()*n);
      for ( int i = 0; i < basis_.localsize(); i++ )
      {
        sum += fac * occ_[nn] * kpg2[i] * norm(p[i]);
      }
    }
  }
  
  ctxt_.dsum(1,1,&sum,1);
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::compute_density(FourierTransform& ft, 
  double weight, double* rho) const
{
  // compute density of the states residing on my column of ctxt_
  assert(occ_.size() == c_.n());
  tmap["density"].start();
  vector<complex<double> > tmp(ft.np012loc());
  
  assert(basis_.cell().volume() > 0.0);
  const double omega_inv = 1.0 / basis_.cell().volume();
  const int np012loc = ft.np012loc();
  
  for ( int i = 0; i < np012loc; i++ )
    rho[i] = 0.0;

  if ( basis_.real() )
  {
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ )
    {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = omega_inv * occ_[nn];
      const double fac2 = omega_inv * occ_[nn+1];
      
      if ( fac1 + fac2 > 0.0 )
      {
        ft.backward(c_.cvalptr(n*c_.mloc()),
                   c_.cvalptr((n+1)*c_.mloc()),&tmp[0]);
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        for ( int i = 0; i < np012loc; i++ )
        {
          const double psi1 = psi[ii];
          const double psi2 = psi[ii+1];
          rho[i] += fac1 * psi1 * psi1 + fac2 * psi2 * psi2;
          ii++; ii++;
        }
      }
    }
    if ( nstloc() % 2 != 0 )
    {
      const int n = nstloc()-1;
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = omega_inv * occ_[nn];
      
      if ( fac1 > 0.0 )
      {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        for ( int i = 0; i < np012loc; i++ )
        {
          const double psi1 = psi[ii];
          rho[i] += fac1 * psi1 * psi1;
          ii++; ii++;
        }
      }
    }
  }
  else
  {
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ )
    {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac = omega_inv * occ_[nn];
      
      if ( fac > 0.0 )
      {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);
        for ( int i = 0; i < np012loc; i++ )
          rho[i] += fac * norm(tmp[i]);
      }
    }
  }
  
  tmap["density"].stop();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::rs_mul_add(FourierTransform& ft, 
  double* v, SlaterDet& sdp) const
{
  // transform states to real space, multiply states by v[r] in real space
  // transform back to reciprocal space and add to sdp
  // sdp[n] += v * sd[n]
  tmap["rs_mul_add"].start();
  
  vector<complex<double> > tmp(ft.np012loc());
  vector<complex<double> > ctmp(2*c_.mloc());
  
  const int np012loc = ft.np012loc();
  const int mloc = c_.mloc();
  double* p = (double*) &tmp[0];
  double* dcp = (double*) sdp.c().valptr();

  if ( basis_.real() )
  {
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ )
    {
      ft.backward(c_.cvalptr(n*mloc),
                 c_.cvalptr((n+1)*mloc),&tmp[0]);
      int ii = 0;
      for ( int i = 0; i < np012loc; i++ )
      {
        const double psi1 = p[ii];
        const double psi2 = p[ii+1];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = vii * psi2;
        ii++; ii++;
      }
      ft.forward(&tmp[0], &ctmp[0], &ctmp[mloc]);
      int len = 4 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy_(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
    if ( nstloc() % 2 != 0 )
    {
      const int n = nstloc()-1;
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);
      int ii = 0;
      for ( int i = 0; i < np012loc; i++ )
      {
        const double psi1 = p[ii];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = 0.0;
        ii++; ii++;
      }
      ft.forward(&tmp[0], &ctmp[0]);
      int len = 2 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy_(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
  }
  else
  {
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ )
    {
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);
      for ( int i = 0; i < np012loc; i++ )
        tmp[i] *= v[i];
      ft.forward(&tmp[0], &ctmp[0]);
      int len = 2 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy_(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
  }
  
  tmap["rs_mul_add"].stop();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::gram(void)
{
  tmap["gram"].start();
  if ( basis_.real() )
  {
    // k = 0 case
    // create a DoubleMatrix proxy for c_
    DoubleMatrix c_proxy(c_);
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.syrk('l','t',2.0,c_proxy,0.0); 
    s.syr('l',-1.0,c_proxy,0,'r');
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
    // solve triangular system X * L^T = C
    c_proxy.trsm('r','l','t','n',1.0,s);
  }
  else
  {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.herk('l','c',1.0,c_,1.0);
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
    // solve triangular system X * L^T = C
    c_.trsm('r','l','c','n',1.0,s);
  }
  tmap["gram"].stop();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::riccati(SlaterDet& sd)
{
  tmap["riccati"].start();
  if ( basis_.real() )
  {
    // k = 0 case
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();
 
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // DoubleMatrix proxy for c_ and sd.c()
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());
    
    // Factor -1.0 in next line: -0.5 from definition of s, 2.0 for G and -G
    s.syrk('l','t',-1.0,c_proxy,0.5); // s = 0.5 * ( I - A )
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',0.5,c_proxy,0,'r');
    // factor -2.0 in next line: G and -G
    r.gemm('t','n',-2.0,sdc_proxy,c_proxy,1.0); // r = ( I - B )
    // rank-1 update using first row of sdc_proxy() and c_proxy
    r.ger(1.0,sdc_proxy,0,c_proxy,0);
 
    xm = s;
    xm.symmetrize('l');
 
    s.syrk('l','t',0.5,r,1.0); // s = s + 0.5 * r^T * r
    s.symmetrize('l');
 
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon )
    {
      // x = s - 0.5 * ( r - xm )^T * ( r - xm )
      // Note: t and r are not symmetric, x, xm, and s are symmetric

      for ( int i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];
 
      x = s;
      x.syrk('l','t',-0.5,t,1.0);
 
      // get full matrix x
      x.symmetrize('l');
 
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];
 
      diff = t.nrm2();
 
      xm = x;
      iter++;
    }
    c_proxy.symm('r','l',1.0,x,sdc_proxy,1.0);
  }
  else
  {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();
 
    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // s = 0.5 * ( I - A )
    s.herk('l','c',-0.5,c_,0.5);
    // r = ( I - B )
    r.gemm('c','n',-1.0,sd.c(),c_,1.0);
 
    xm = s;
    xm.symmetrize('l');
 
    // s = s + 0.5 * r^H * r
    s.herk('l','c',0.5,r,1.0);
    s.symmetrize('l');
 
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;
 
    while ( iter < maxiter && diff > epsilon )
    {
      // x = s - 0.5 * ( r - xm )^H * ( r - xm )
      // Note: t and r are not hermitian, x, xm, and s are hermitian

      for ( int i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];
 
      x = s;
      x.herk('l','c',-0.5,t,1.0);
      x.symmetrize('l');
 
      for ( int i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];
 
      diff = t.nrm2();
 
      xm = x;
      iter++;
    }
    c_.hemm('r','l',1.0,x,sd.c(),1.0);
  }
  tmap["riccati"].stop();
}
  
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(int nel, int nspin)
{
  // compute occupation numbers as 0.0, 1.0 or 2.0
  // if nspin = 1: use 0, 1 or 2
  // if nspin = 2: use 0 or 1;
  assert (nel >= 0);
  assert (occ_.size() == c_.n());
  if ( nspin == 1 )
  {
    assert (nel <= 2*c_.n());
    int ndouble = nel/2;
    for ( int n = 0; n < ndouble; n++ )
      occ_[n] = 2.0;
    for ( int n = ndouble; n < ndouble+nel%2; n++ )
      occ_[n] = 1.0;
    for ( int n = ndouble+nel%2; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else if ( nspin == 2 )
  {
    assert (nel <= c_.n());
    for ( int n = 0; n < nel; n++ )
      occ_[n] = 1.0;
    for ( int n = nel; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::total_charge(void)
{
  // compute total charge from occ_[i]
  double sum = 0.0;
  for ( int n = 0; n < occ_.size(); n++ )
  {
    sum += occ_[n];
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(int nspin, double mu, double temp)
{
  // compute occupation numbers using a Fermi distribution f(mu,temp)
  // and the eigenvalues in eig_[i]
  assert(nspin==1 || nspin==2);
  assert (occ_.size() == c_.n());
  assert (eig_.size() == c_.n());
  if ( nspin == 1 )
  {
    for ( int n = 0; n < eig_.size(); n++ )
    {
      occ_[n] = 2.0 * fermi(eig_[n],mu,temp);
    }
  }
  else if ( nspin == 2 )
  {
    for ( int n = 0; n < eig_.size(); n++ )
    {
      occ_[n] = fermi(eig_[n],mu,temp);
    }
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::fermi(double e, double mu, double fermitemp)
{
  // e, mu in Hartree, fermitemp in Kelvin

  if ( fermitemp == 0.0 )
  {
    if ( e < mu ) return 1.0;
    else if ( e == mu ) return 0.5;
    else return 0.0;
  }
  const double kb = 3.1667907e-6; // Hartree/Kelvin
  const double kt = kb * fermitemp;
  double arg = ( e - mu ) / kt;

  if ( arg < -30.0 ) return 1.0;
  if ( arg >  30.0 ) return 0.0;

  return 1.0 / ( 1.0 + exp ( arg ) );
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::entropy(int nspin)
{
  // return dimensionless entropy
  // the contribution to the free energy is - t_kelvin * k_boltz * wf.entropy()

  assert(nspin==1 || nspin==2);
  const double fac = ( nspin > 1 ) ? 1.0 : 2.0;
  double sum = 0.0;
  for ( int n = 0; n < occ_.size(); n++ )
  {
    const double f = occ_[n] / fac;
    if ( f > 0.0  &&  f < 1.0 )
    {
      sum -= fac * ( f * log(f) + (1.0-f) * log(1.0-f) );
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::ortho_error(void)
{
  // deviation from orthogonality of c_
  tmap["ortho_error"].start();
  double error;
  if ( basis_.real() )
  {
    // k = 0 case
    // declare a proxy DoubleMatrix for c_
    DoubleMatrix c_proxy(c_);
    
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    
    // real symmetric rank-k update
    // factor 2.0 in next line: G and -G
    s.syrk('l','t',2.0,c_proxy,0.0); // compute real overlap matrix
 
    // correct for double counting of G=0
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',-1.0,c_proxy,0,'r');
 
    DoubleMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();
    
    s -= id; // subtract identity matrix from S
    
    error = s.nrm2();
  }
  else
  {
    tmap["ortho_error"].start();
    // k != 0 case
    
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.herk('l','c',1.0,c_,0.0);
    
    ComplexMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();
 
    s -= id; // subtract identity matrix from S
    
    error = s.nrm2();
  }
  tmap["ortho_error"].stop();
  return error;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::randomize(double amplitude)
{
  // Note: randomization results depend on the process grid size and shape
  srand48(ctxt_.myproc());
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    for ( int i = 0; i < basis_.localsize(); i++ )
    {
      double re = drand48();
      double im = drand48();
      p[i] += amplitude * complex<double>(re,im);
    }
  }
  cleanup();
  gram();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::cleanup(void)
{
  // set Im( c(G=0) ) to zero and 
  // set the empty rows of the matrix c_ to zero
  // The empty rows are located between i = basis_.localsize() and 
  // c_.mloc(). Empty rows are necessary to insure that the 
  // local size c_.mloc() is the same on all processes, while the 
  // local basis size is not.
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    // reset imaginary part of G=0 component to zero
    if ( ctxt_.myrow() == 0 )
    {
      // index of G=0 element
      int izero;
      if ( basis_.real() )
        izero = 0;
      else
        izero = basis_.rod_size(0)/2;
      //cout << " izero = " << izero << " G = " << basis_.kv(3*izero) << " " 
      //     << basis_.kv(3*izero+1) << " " << basis_.kv(3*izero+2) << endl;
      p[izero] = complex<double> ( p[izero].real(), 0.0);
    }
    // reset values of empty rows of c_ to zero
    for ( int i = basis_.localsize(); i < c_.mloc(); i++ )
    {
      p[i] = 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
SlaterDet& SlaterDet::operator=(SlaterDet& rhs)
{
  if ( this == &rhs ) return *this;
  assert(ctxt_.ictxt() == rhs.context().ictxt());
  c_ = rhs.c_;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::memsize(void) const
{
  return basis_.memsize() + c_.memsize();
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::localmemsize(void) const
{
  return basis_.localmemsize() + c_.localmemsize();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::print(ostream& os, string encoding)
{
  FourierTransform ft(basis_,basis_.np(0),basis_.np(1),basis_.np(2));
  vector<complex<double> > wftmp(ft.np012loc());
  vector<double> wftmpr(ft.np012());
  
  if ( ctxt_.onpe0() )
  {
    const double weight = 1.0; //!! fixed determinant weight to 1.0
    //!! no spin attribute written
    os << "<slater_determinant kpoint=\"" << basis_.kpoint() << "\"\n"
       << "  weight=\"" << weight << "\""
       << " size=\"" << nst() << "\">" << endl;
 
    os << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">" 
       << endl;
    os.setf(ios::fixed,ios::floatfield);
    os.setf(ios::right,ios::adjustfield);
    for ( int i = 0; i < nst(); i++ )
    {
      os << " " << setprecision(8) << occ_[i];
      if ( i%10 == 9 )
        os << endl;
    }
    if ( nst()%10 != 0 )
      os << endl;
    os << "</density_matrix>" << endl;
  }
  
  for ( int n = 0; n < nst(); n++ )
  {
    // check if state n resides on mype
    if ( c_.pc(n) == ctxt_.mycol() )
    {
      //cout << " state " << n << " is stored on column " 
      //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
      int nloc = c_.y(n); // local index
      ft.backward(c_.cvalptr(c_.mloc()*nloc),&wftmp[0]);
      
      double *a = (double*) &wftmp[0];
      for ( int i = 0; i < ft.np012loc(); i++ )
        wftmpr[i] = a[2*i];
        
      for ( int i = 0; i < ctxt_.nprow(); i++ )
      {
        if ( i == ctxt_.myrow() )
        {
          int size = ft.np012loc();
          //cout << " process " << ctxt_.mype() << " sending block " << i
          //     << " of state "
          //     << n << " to task 0, size = " << size << endl;
          ctxt_.isend(1,1,&size,1,0,0);
          ctxt_.dsend(size,1,&wftmpr[0],1,0,0);
        }
      }
    }
    if ( ctxt_.onpe0() )
    {
      for ( int i = 0; i < ctxt_.nprow(); i++ )
      {
        int size = 0;
        ctxt_.irecv(1,1,&size,1,i,c_.pc(n));
        int istart = ft.np0() * ft.np1() * ft.np2_first(i);
        //cout << " task 0 (" << ctxt_.mype() << ") receiving block " << i
        //     << " of state "
        //     << n << " size=" << size << " istart=" << istart << endl;
        ctxt_.drecv(size,1,&wftmpr[istart],1,i,c_.pc(n));
      }
      
      // wftmpr is now complete on task 0

      if ( encoding == "base64" )
      {
        #if AIX
        byteswap_double(ft.np012(),&wftmpr[0]);
        #endif
        unsigned int outlen;
        XMLByte* b = Base64::encode((XMLByte*) &wftmpr[0],
          ft.np012()*sizeof(double), &outlen);
        assert(b!=0);
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"base64\">" << endl;
        os.write((char*) b, outlen);
        os << "</grid_function>\n";
        XMLString::release(&b);
      }
      else
      {
        // encoding == "text" or unknown encoding
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"text\">" << endl;
        int count = 0;
        for ( int k = 0; k < ft.np2(); k++ )
          for ( int j = 0; j < ft.np1(); j++ )
            for ( int i = 0; i < ft.np0(); i++ )
            {
              os << " " << wftmpr[ft.index(i,j,k)];
              if ( count++%4 == 3)
                os << "\n";
            }
        if ( count%4 != 0 )
          os << "\n";
        os << "</grid_function>\n";
      }

    }
  }
  if ( ctxt_.onpe0() )
    os << "</slater_determinant>" << endl;
}

#if USE_CSTDIO_LFS
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::write(FILE* outfile, string encoding)
{
  FourierTransform ft(basis_,basis_.np(0),basis_.np(1),basis_.np(2));
  vector<complex<double> > wftmp(ft.np012loc());
  vector<double> wftmpr(ft.np012());
  
  ostringstream os;
  
  if ( ctxt_.onpe0() )
  {
    const double weight = 1.0; //!! fixed determinant weight to 1.0
    //!! no spin attribute written
    os << "<slater_determinant kpoint=\"" << basis_.kpoint() << "\"\n"
       << "  weight=\"" << weight << "\""
       << " size=\"" << nst() << "\">" << endl;
 
    os << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">" 
       << endl;
    os.setf(ios::fixed,ios::floatfield);
    os.setf(ios::right,ios::adjustfield);
    for ( int i = 0; i < nst(); i++ )
    {
      os << " " << setprecision(8) << occ_[i];
      if ( i%10 == 9 )
        os << endl;
    }
    if ( nst()%10 != 0 )
      os << endl;
    os << "</density_matrix>" << endl;
    
    string str(os.str());
    off_t len = str.length();
    fwrite(str.c_str(),sizeof(char),len,outfile);
  }
  
  for ( int n = 0; n < nst(); n++ )
  {
    // check if state n resides on mype
    if ( c_.pc(n) == ctxt_.mycol() )
    {
      //cout << " state " << n << " is stored on column " 
      //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
      int nloc = c_.y(n); // local index
      ft.backward(c_.cvalptr(c_.mloc()*nloc),&wftmp[0]);
      
      double *a = (double*) &wftmp[0];
      for ( int i = 0; i < ft.np012loc(); i++ )
        wftmpr[i] = a[2*i];
        
      for ( int i = 0; i < ctxt_.nprow(); i++ )
      {
        if ( i == ctxt_.myrow() )
        {
          int size = ft.np012loc();
          //cout << " process " << ctxt_.mype() << " sending block " << i
          //     << " of state "
          //     << n << " to task 0, size = " << size << endl;
          ctxt_.isend(1,1,&size,1,0,0);
          ctxt_.dsend(size,1,&wftmpr[0],1,0,0);
        }
      }
    }
    if ( ctxt_.onpe0() )
    {
      for ( int i = 0; i < ctxt_.nprow(); i++ )
      {
        int size = 0;
        ctxt_.irecv(1,1,&size,1,i,c_.pc(n));
        int istart = ft.np0() * ft.np1() * ft.np2_first(i);
        //cout << " task 0 (" << ctxt_.mype() << ") receiving block " << i
        //     << " of state "
        //     << n << " size=" << size << " istart=" << istart << endl;
        ctxt_.drecv(size,1,&wftmpr[istart],1,i,c_.pc(n));
      }
      
      // wftmpr is now complete on task 0

      if ( encoding == "base64" )
      {
        #if AIX
        byteswap_double(ft.np012(),&wftmpr[0]);
        #endif
        unsigned int outlen;
        XMLByte* b = Base64::encode((XMLByte*) &wftmpr[0],
          ft.np012()*sizeof(double), &outlen);
        assert(b!=0);
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os.str("");
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"base64\">" << endl;
        string str(os.str());
        off_t len = str.length();
        fwrite(str.c_str(),sizeof(char),len,outfile);
        fwrite((char*) b,sizeof(char),outlen,outfile);
        //os.write((char*) b, outlen);
        os.str("");
        os << "</grid_function>\n";
        str = os.str();
        len = str.length();
        fwrite(str.c_str(),sizeof(char),len,outfile);
        XMLString::release(&b);
      }
      else
      {
        // encoding == "text" or unknown encoding
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os.str("");
        os << "<grid_function type=\"double\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"text\">" << endl;
        string str(os.str());
        off_t len = str.length();
        fwrite(str.c_str(),sizeof(char),len,outfile);
        int count = 0;
        for ( int k = 0; k < ft.np2(); k++ )
        {
          os.str("");
          for ( int j = 0; j < ft.np1(); j++ )
          {
            for ( int i = 0; i < ft.np0(); i++ )
            {
              os << " " << wftmpr[ft.index(i,j,k)];
              if ( count++%4 == 3)
                os << "\n";
            }
          }
          str = os.str();
          len = str.length();
          fwrite(str.c_str(),sizeof(char),len,outfile);
        }
        os.str("");
        if ( count%4 != 0 )
          os << "\n";
        os << "</grid_function>\n";
        str = os.str();
        len = str.length();
        fwrite(str.c_str(),sizeof(char),len,outfile);
      }

    }
  }
  
  if ( ctxt_.onpe0() )
  {
    os.str("");
    os << "</slater_determinant>" << endl;
    string str(os.str());
    off_t len = str.length();
    fwrite(str.c_str(),sizeof(char),len,outfile);
  }
}
#endif

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::info(ostream& os)
{  
  if ( ctxt_.onpe0() )
  {
    const double weight = 1.0; //!! fixed determinant weight to 1.0
    //!! no spin attribute written
    os << "<slater_determinant kpoint=\"" << basis_.kpoint() << "\"\n"
       << "  weight=\"" << weight << "\""
       << " size=\"" << nst() << "\">" << endl;
    os << " <!-- basis size: " << basis_.size() << " -->" << endl;
    os << "  <!-- c dimensions: "
       << c_.m() << "x" << c_.n()
       << "   (" << c_.mb() << "x" << c_.nb() << " blocks)" << " -->" << endl;
    os << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">" 
       << endl;
    os << "</density_matrix>" << endl;
    os << "</slater_determinant>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, SlaterDet& sd)
{
  sd.print(os,"text");
  return os;
}
