////////////////////////////////////////////////////////////////////////////////
//
// testLineMinimizer.C
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "LineMinimizer.h"
using namespace std;

const double fac = 0.1;
#if 1
double ff(double x) { return fac*(x*x*x*x -100*x*x - 100*x); }
double gg(double x) { return fac*(4*x*x*x - 200*x - 100); }
#else
double ff(double x) { return fac*(x*x*x*x); }
double gg(double x) { return fac*(4*x*x*x); }
#endif
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // use: testLineMinimizer niter xstart sigma1 sigma2 alpha_start
  if ( argc != 6 )
  {
    cout << "use: testLineMinimizer niter xstart sigma1 sigma2 alpha_start"
         << endl;
    return 1;
  }

  double x,x0;
  const int n = atoi(argv[1]);
  x = atof(argv[2]);
  const double sigma1 = atof(argv[3]);
  const double sigma2 = atof(argv[4]);
  const double alpha_start = atof(argv[5]);
  x0 = x;
  double alpha = 0;
  LineMinimizer linmin;
  linmin.set_sigma1(sigma1);
  linmin.set_sigma2(sigma2);
  linmin.set_alpha_start(alpha_start);

  // descent direction p
  double p = -gg(x0);

  bool done = false;
  int i = 0;
  while ( !done && i < n )
  {
    i++;
    double f = ff(x);
    double g = gg(x);
    // fp = derivative in the descent direction p
    double fp = g * p;
    done = linmin.done();
    cout << "i=" << i << " x=" << x << " f=" << f << " fp=" << fp << endl;
    if ( !done )
    {
      alpha = linmin.next_alpha(alpha,f,fp);
      cout << "alpha=" << alpha << endl;
      x = x0 + alpha * p;
    }
  }
}
