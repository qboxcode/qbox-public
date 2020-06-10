////////////////////////////////////////////////////////////////////////////////
//
// testVWN.cpp
//
// test the Vosko-Wilk-Nusair functional
//
////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include "VWNFunctional.h"
#include "LDAFunctional.h"
#include <cmath>
using namespace std;

// c1 = (3.D0/(4.D0*pi))**third
const double c1 = 0.6203504908994001;
const double c3 = -0.610887057711;

double rho_of_rs(double rs)
{
  return pow(c1/rs,3.0);
}

int main()
{
  const int n = 500;
  vector<vector<double> > rho;
  rho.resize(1);
  rho[0].resize(n);

  LDAFunctional f(rho);

  double rs[] = { 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0,
                  15.0, 20.0, 50.0, 100.0 };

  // print table for given rs values (Table 5 in paper)
  for ( int i = 0; i < 13; i++ )
    rho[0][i] = rho_of_rs(rs[i]);

  f.setxc();
  for ( int i = 0; i < 13; i++ )
  {
    double vx = c3 / rs[i];
    double ex = 0.75 * vx;
    cout << "rs=" << rs[i] << " " << f.rho[i] << " "
         << f.exc[i] << " " << f.vxc1[i] << endl;
  }
  cout << endl;

#if 0
  const double drh = 0.01;
  for ( int i = 0; i < rho[0].size(); i++ )
  {
    double rh = drh * i;
    rho[0][i] = rh;
  }
  f.setxc();
  for ( int i = 1; i < rho[0].size()-1; i++ )
  {
    cout << f.rho[i] << " " << f.exc[i] << " " << f.vxc1[i]
         << " " << (f.rho[i+1]*f.exc[i+1]-f.rho[i-1]*f.exc[i-1])/(2*drh)
         << endl;
  }
#endif

  // spin-polarized test
  rho.resize(2);
  rho[0].resize(500);
  rho[1].resize(500);
  for ( int i = 0; i < 13; i++ )
  {
    rho[0][i] = rho_of_rs(rs[i]);
    rho[1][i] = 0.0;
  }

  LDAFunctional fs(rho);

  fs.setxc();
  for ( int i = 0; i < 13; i++ )
  {
    const double c4 = -0.769669463118;
    double vx = c4 / rs[i];
    double ex = 0.75 * vx;
    cout << "rs=" << rs[i] << " " << fs.rho_up[i] << " "
         << fs.exc[i] << " " << fs.vxc1_up[i] << endl;
  }
  cout << endl;

}
