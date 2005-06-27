////////////////////////////////////////////////////////////////////////////////
//
// ConstraintSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintSet.h,v 1.1 2005-06-27 22:34:46 fgygi Exp $

#ifndef CONSTRAINTSET_H
#define CONSTRAINTSET_H

#include <vector>
#include <string>
using namespace std;

class Atom;
class AtomSet;
class Constraint;
class Context;

class ConstraintSet
{
  private:
  
  const Context& ctxt_;
  vector<Constraint *> constraint_list;
  
  public:

  ConstraintSet(const Context& ctxt) : ctxt_(ctxt) {}
  bool set_constraint(AtomSet &atoms, int argc, char **argv);
  bool del_constraint(int argc, char **argv);
  void list_constraints(ostream &os);
  int size(void) const { return constraint_list.size(); }
  void enforce(AtomSet& atoms);
  void enforce_r(const vector<vector<double> > &r0,
                 vector<vector<double> > &rp);
  void enforce_v(const vector<vector<double> > &r0,
                 vector<vector<double> > &v0);
  double projection(const vector<vector<double> > &r0,
                    const vector<vector<double> > &x);
  void update_constraints(double dt);
  void setup(AtomSet& atoms);
};
#endif
