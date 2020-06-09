////////////////////////////////////////////////////////////////////////////////
//
// testCGOptimizer.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include "CGOptimizer.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Rosenbrock function
double rb(const valarray<double> &x, valarray<double> &df)
{
  int n = x.size();
  assert ( n%2 == 0 );
  double sum = 0.0;
  for ( int i = 0; i < n; i++, i++ )
  {
    sum += ( ( 1.0 - x[i] )*( 1.0 - x[i] ) +
             100.0 * ( x[i+1] - x[i]*x[i] ) * ( x[i+1] - x[i]*x[i] ) );
    df[i+1] = 200.0 * ( x[i+1] - x[i]*x[i] );
    df[i] = -2.0 * ( x[i]*df[i+1] + ( 1.0 - x[i] ) );
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double quad1(const valarray<double> &x, valarray<double> &df)
{
  df[0] = 2 * x[0];
  return x[0]*x[0];
}

////////////////////////////////////////////////////////////////////////////////
double quad2(const valarray<double> &x, valarray<double> &df)
{
  df[0] = 2 * x[0];
  df[1] = 4 * x[1];
  return x[0]*x[0] + 2.0*x[1]*x[1];
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if ( argc != 3 )
  {
    cout << "use: testCGOptimizer ndim niter" << endl;
    return 1;
  }
  const int n = atoi(argv[1]);
  const int niter = atoi(argv[2]);

  CGOptimizer cgop(n);
  cgop.set_alpha_start(0.01);
  cgop.set_beta_max(20);
  cgop.set_debug_print();

  valarray<double> x(n), xp(n), g(n);
  double f;

  for ( int i = 0; i < n; i++ )
    x[i] = 0.01*(i+1);

  for ( int i = 0; i < niter; i++ )
  {
    f = rb(x,g);
    //f = quad2(x,g);
    cout << "iteration " << i << " f=" << f << endl;
    for ( int i = 0; i < n; i++ )
      cout << x[i] << " " << g[i] << endl;

    cgop.compute_xp(x,f,g,xp);
    x = xp;
  }
}
