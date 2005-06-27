////////////////////////////////////////////////////////////////////////////////
//
//  Constraint.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Constraint.h,v 1.1 2005-06-27 22:34:46 fgygi Exp $

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <string>
#include <vector>
#include <cassert>
using namespace std;

class AtomSet;

class Constraint
{
  protected:
  
  vector<string> names_; // names of atoms involved in the constraint
  
  public:
  
  virtual string type(void) const = 0;
  virtual double value(void) const = 0;
  virtual double velocity(void) const = 0;
  virtual double tolerance(void) const = 0;
  virtual void set_value(double value) = 0;
  virtual void set_velocity(double velocity) = 0;
  virtual bool enforce_r(const vector<vector<double> > &r0,
                         vector<vector<double> > &rp) const = 0;
  virtual bool enforce_v(const vector<vector<double> > &r0,
                         vector<vector<double> > &v0) const = 0;
  virtual double projection(const vector<vector<double> > &r0,
                            const vector<vector<double> > &x) const = 0;
  virtual void update(double dt) = 0;
  virtual void setup(const AtomSet& atoms) = 0;
  virtual ostream& print(ostream &os) = 0;
  string names(int i)
  {
    assert( i >= 0 && i < names_.size() );
    return names_[i];
  }
};
ostream& operator << ( ostream &os, Constraint &c );
#endif
  
  
