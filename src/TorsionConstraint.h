////////////////////////////////////////////////////////////////////////////////
//
//  TorsionConstraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: TorsionConstraint.h,v 1.1 2005-06-27 22:34:46 fgygi Exp $

#ifndef TORSIONCONSTRAINT_H
#define TORSIONCONSTRAINT_H

#include "Constraint.h"
#include "D3vector.h"
#include <cassert>
class AtomSet;

class TorsionConstraint : public Constraint
{
  string name1_, name2_, name3_, name4_;
  int    ia1_, ia2_, ia3_, ia4_, is1_, is2_, is3_, is4_;
  double m1_, m2_, m3_, m4_, m1_inv_, m2_inv_, m3_inv_, m4_inv_;
  double angle_, velocity_, tol_, sin_angle_, cos_angle_;
  double sigma(D3vector a, D3vector b,
               D3vector c, D3vector d) const;
  void grad_sigma(const D3vector &r1, const D3vector &r2,
                    const D3vector &r3, const D3vector &r4,
                    D3vector &g1, D3vector &g2,D3vector &g3,D3vector &g4) const;
  double torsion_angle(D3vector a, D3vector b,
                       D3vector c, D3vector d) const;
  
  public:
  
  TorsionConstraint(string name1, string name2, string name3, string name4,
                     double angle, double velocity, double tolerance):
  name1_(name1), name2_(name2), name3_(name3), name4_(name4),
  velocity_(velocity),
  tol_(tolerance), m1_(0.0), m2_(0.0), m3_(0.0), m4_(0.0)
  {
    set_value(angle);
    names_.resize(4);
    names_[0] = name1_;
    names_[1] = name2_;
    names_[2] = name3_;
    names_[3] = name4_;
  }
  
  string type(void) const { return "torsion"; }
  double value(void) const { return angle_; }
  double velocity(void) const { return velocity_; }
  double tolerance(void) const { return tol_; }
  void set_value(double value)
  {
    angle_ = value;
    if ( angle_ < -180.0 ) angle_ = 180.0;
    if ( angle_ >  180.0 ) angle_ = 180.0;
    sin_angle_ = sin((M_PI/180.0)*angle_);
    cos_angle_ = cos((M_PI/180.0)*angle_);
  }
  void set_velocity(double velocity)
  {
    velocity_ = velocity;
  }
  
  void setup(const AtomSet& atoms);
  void update(double dt);
  bool enforce_r(const vector<vector<double> > &r0,
                 vector<vector<double> > &rp) const;
  bool enforce_v(const vector<vector<double> > &r0,
                 vector<vector<double> > &v0) const;
  double projection(const vector<vector<double> > &r0,
                    const vector<vector<double> > &x) const;
  ostream& print( ostream& os );
  
};
#endif
