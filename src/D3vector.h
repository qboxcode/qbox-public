////////////////////////////////////////////////////////////////////////////////
//
// D3vector.h
//
// double 3-vectors
//
////////////////////////////////////////////////////////////////////////////////
// $ Id: $

#ifndef D3VECTOR_H
#define D3VECTOR_H
#include <iostream>
#include <cmath>
#include <cassert>
using namespace std;

class D3vector
{
  public:

  double x, y, z; 

  // explicit constructor to avoid implicit conversion from double to D3vector
  explicit D3vector(const double& xv = 0, const double& yv = 0, const double& zv = 0) :
  x(xv), y(yv), z(zv) {}
  
  double& operator[](const int &i)
  {
    assert(i>=0 && i <3);
    if ( i == 0 ) return x;
    else if ( i == 1 ) return y;
    else return z;
  }
  
  bool operator==(const D3vector &rhs) const
  {
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
  
  bool operator!=(const D3vector &rhs) const
  {
    return x != rhs.x || y != rhs.y || z != rhs.z;
  }

  D3vector& operator += ( const D3vector& rhs )
  {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }

  D3vector& operator -= ( const D3vector& rhs )
  {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *this;
  }

  D3vector& operator *= ( const double& rhs )
  {
    x *= rhs; y *= rhs; z *= rhs;
    return *this;
  }

  D3vector& operator /= ( const double& rhs )
  {
    x /= rhs; y /= rhs; z /= rhs;
    return *this;
  }

  friend const D3vector operator + (const D3vector& lhs, const D3vector& rhs )
  {
    return D3vector(lhs) += rhs;
  }

  friend const D3vector operator - ( const D3vector& a, const D3vector& b )
  {
    return D3vector(a) -= b;
  }

  friend D3vector operator - ( const D3vector& a ) // unary minus
  {
    return D3vector( -a.x, -a.y, -a.z );
  }

  friend D3vector operator * ( const double& a, const D3vector& b )
  {
    return D3vector(b) *= a;
  }

  friend D3vector operator * ( const D3vector& a, const double& b )
  {
    return D3vector(a) *= b;
  }

  friend D3vector operator / ( const D3vector& a, const double& b )
  {
    return D3vector(a) /= b;
  }

  // scalar product
  friend double operator * ( const D3vector& a, const D3vector& b )
  {
    return a.x * b.x + a.y * b.y + a.z * b.z ;
  }

  friend D3vector operator ^ ( const D3vector& a, const D3vector& b )
  {
    return D3vector( a.y * b.z - a.z * b.y ,
                     a.z * b.x - a.x * b.z ,
                     a.x * b.y - a.y * b.x );
  }

  friend D3vector rotate ( const D3vector& x, const D3vector& w )
  {
    if ( length(x) == 0.0 ) return x; // x has zero length
    double theta = length( w );       // rotate by zero
    if ( theta == 0.0 ) return x;
    D3vector ew = normalized ( w );
    D3vector v = w ^ x;
    if ( length( v ) == 0.0 ) return x; // x is parallel to the rotation axis
    v = normalized( v );
    D3vector u = v ^ ew;
    double p = x * u;
    return  (x*ew)*ew + p*cos(theta)*u + p*sin(theta)*v ;
  }

  friend double length( const D3vector& a )
  {
    return sqrt( a.x * a.x + a.y * a.y + a.z * a.z );
  }

  friend double norm( const D3vector& a )
  {
    return a.x * a.x + a.y * a.y + a.z * a.z;
  }

  friend D3vector normalized( const D3vector a )
  {
    return a / length( a );
  }

  friend ostream& operator << ( ostream& os, const D3vector& v )
  {
    os << v.x << " " << v.y << " " << v.z;
    return os;
  }

  friend istream& operator >> ( istream& is, D3vector& v )
  {
    is >> v.x >> v.y >> v.z ;
    return is;
  }

};
#endif
