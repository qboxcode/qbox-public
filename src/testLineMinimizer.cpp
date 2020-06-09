////////////////////////////////////////////////////////////////////////////////
//
// testLineMinimizer.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cassert>
#include <cstdlib>
#include "LineMinimizer.h"
using namespace std;

const double fac = 0.1;
void ff(int ifun, double x, double *f, double *g)
{
  const double s = 0.45;
  switch (ifun)
  {
    case 1:
      *f = fac*(x*x*x*x -100*x*x - 100*x);
      *g = fac*(4*x*x*x - 200*x - 100);
      break;

    case 2:
      *f = fac*(x*x*x*x);
      *g = fac*(4*x*x*x);
      break;

    case 3:
      // function (3.1) of More and Thuente
      if ( x <= 1.0 )
      {
        double s = 0.3;
        *f = 0.5*(1.0-s)*x*x - x;
        *g = (1.0-s)*x - 1.0;
      }
      else
      {
        *f = 0.5*(s-1.0) - s*x;
        *g = -s;
      }
      break;

    case 4:
      // function (3.3) of More and Thuente
      *f = -3*x/(x*x+2.0) - 0.03*x;
      *g = (-3*(x*x+2.0) + 6*x*x)/((x*x+2)*(x*x+2)) - 0.03;
      break;

    default:
      cout << "incorrect ifun value" << endl;
      assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // use: testLineMinimizer ifun niter xstart sigma1 sigma2 alpha_start
  if ( argc != 7 )
  {
    cout << "use:testLineMinimizer ifun, niter xstart sigma1 sigma2 alpha_start"
         << endl;
    return 1;
  }

  const int ifun = atoi(argv[1]);
  const int n = atoi(argv[2]);
  double x = atof(argv[3]);
  const double sigma1 = atof(argv[4]);
  const double sigma2 = atof(argv[5]);
  const double alpha_start = atof(argv[6]);
  double x0 = x;
  double f,g;
  double alpha = 0;
  LineMinimizer linmin;
  linmin.set_sigma1(sigma1);
  linmin.set_sigma2(sigma2);
  linmin.set_alpha_start(alpha_start);
  linmin.set_debug_print();

  // descent direction p
  ff(ifun,x0,&f,&g);
  double p = -g;

  bool done = false;
  int i = 0;
  while ( !done && i < n )
  {
    cout << "========================================================" << endl;
    i++;
    ff(ifun,x,&f,&g);
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
  cout << "========================================================" << endl;
}
