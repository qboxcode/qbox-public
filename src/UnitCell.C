////////////////////////////////////////////////////////////////////////////////
//
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: UnitCell.C,v 1.7 2004-02-04 19:55:17 fgygi Exp $

#include "UnitCell.h"
#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void UnitCell::set(const D3vector& a0, const D3vector& a1, const D3vector& a2)
{
  a_[0] = a0; a_[1] = a1, a_[2] = a2;
  volume_ = a0 * ( a1 ^ a2 );
  if ( volume_ > 0.0 )
  {
    double fac = 2.0 * M_PI / volume_;
    b_[0] = fac * a1 ^ a2;
    b_[1] = fac * a2 ^ a0;
    b_[2] = fac * a0 ^ a1;
    
  }
  else
  {
    b_[0] = b_[1] = b_[2] = D3vector(0.0,0.0,0.0);
  }
  
  an_[0]  = a_[0];
  an_[1]  = a_[1];
  an_[2]  = a_[2];
  an_[3]  = a_[0] + a_[1];
  an_[4]  = a_[0] - a_[1];
  an_[5]  = a_[1] + a_[2];
  an_[6]  = a_[1] - a_[2];
  an_[7]  = a_[2] + a_[0];
  an_[8]  = a_[2] - a_[0];
  an_[9]  = a_[0] + a_[1] + a_[2];
  an_[10] = a_[0] - a_[1] - a_[2];
  an_[11] = a_[0] + a_[1] - a_[2];
  an_[12] = a_[0] - a_[1] + a_[2];
 
  bn_[0]  = b_[0];
  bn_[1]  = b_[1];
  bn_[2]  = b_[2];
  bn_[3]  = b_[0] + b_[1];
  bn_[4]  = b_[0] - b_[1];
  bn_[5]  = b_[1] + b_[2];
  bn_[6]  = b_[1] - b_[2];
  bn_[7]  = b_[2] + b_[0];
  bn_[8]  = b_[2] - b_[0];
  bn_[9]  = b_[0] + b_[1] + b_[2];
  bn_[10] = b_[0] - b_[1] - b_[2];
  bn_[11] = b_[0] + b_[1] - b_[2];
  bn_[12] = b_[0] - b_[1] + b_[2];
 
  for ( int i = 0; i < 13; i++ )
  {
    an2h_[i] = 0.5 * norm(an_[i]);
    bn2h_[i] = 0.5 * norm(bn_[i]);
  }
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::in_ws(const D3vector& v) const
{
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( abs(v*an_[i]) <= an2h_[i] ) ;
    i++;
  }
  return in;
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::fold_in_ws(D3vector& v) const
{
  bool done = false;
  while ( !done )
  {
    done = true;
    for ( int i = 0; (i < 13) && done; i++ )
    {
      const double sp = v*an_[i];
      if ( sp > an2h_[i] )
      {
        done = false;
        do
          v -= an_[i];
        while ( v*an_[i] > an2h_[i] );
      }
      else if ( sp < -an2h_[i] )
      {
        done = false;
        do
          v += an_[i];
        while ( v*an_[i] < -an2h_[i] );
      }
    }
  }
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::in_bz(const D3vector& k) const
{
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( abs(k*bn_[i]) <= bn2h_[i] ) ;
    i++;
  }
  return in;
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::fold_in_bz(D3vector& k) const
{
  bool done = false;
  while ( !done )
  {
    done = true;
    for ( int i = 0; (i < 13) && done; i++ )
    {
      double sp = k*bn_[i];
      if ( sp > bn2h_[i] )
      {
        done = false;
        do
          k -= bn_[i];
        while ( k*bn_[i] > bn2h_[i] );
      }
      else if ( sp < -bn2h_[i] )
      {
        done = false;
        do
          k += bn_[i];
        while ( k*bn_[i] < -bn2h_[i] );
      }
    }
  }
}
 
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::encloses(const UnitCell& c) const
{
  D3vector center(0.5*an_[9]); // an_[9] is a_[0] + a_[1] + a_[2]
  bool in = true;
  int i = 0;
  while ( i < 13 && in )
  {
    in = ( contains(center+0.5*c.an_[i]) );
    i++;
  }
  return in;
}

////////////////////////////////////////////////////////////////////////////////
bool UnitCell::contains(D3vector v) const
{
  const double fac = 0.5 / ( 2.0 * M_PI );
  const double p0 = fac * v * b_[0];
  const double p1 = fac * v * b_[1];
  const double p2 = fac * v * b_[2];
  return ( (p0 > 0.0) && (p0 <= 1.0) &&
           (p1 > 0.0) && (p1 <= 1.0) &&
           (p2 > 0.0) && (p2 <= 1.0) );
}
 
////////////////////////////////////////////////////////////////////////////////
void UnitCell::print(ostream& os) const
{
  os.setf(ios::fixed,ios::floatfield);
  os << "    <a> " << setw(10) << setprecision(4) << a_[0].x
                   << setw(10) << setprecision(4) << a_[0].y
                   << setw(10) << setprecision(4) << a_[0].z << " </a>" << endl;
  os << "    <b> " << setw(10) << setprecision(4) << a_[1].x
                   << setw(10) << setprecision(4) << a_[1].y
                   << setw(10) << setprecision(4) << a_[1].z << " </b>" << endl;
  os << "    <c> " << setw(10) << setprecision(4) << a_[2].x
                   << setw(10) << setprecision(4) << a_[2].y
                   << setw(10) << setprecision(4) << a_[2].z << " </c>" << endl;
}
  
////////////////////////////////////////////////////////////////////////////////
bool UnitCell::operator==(const UnitCell& c) const
{
  return ( a_[0]==c.a_[0] && a_[1]==c.a_[1] && a_[2]==c.a_[2] );
}
 
////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const UnitCell& cell )
{ 
  cell.print(os); 
  return os;
}
