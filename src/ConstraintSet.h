////////////////////////////////////////////////////////////////////////////////
//
// ConstraintSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ConstraintSet.h,v 1.4 2007-10-19 16:24:04 fgygi Exp $

#ifndef CONSTRAINTSET_H
#define CONSTRAINTSET_H

#include <vector>
#include <string>

class Atom;
class AtomSet;
class Constraint;
class Context;

class ConstraintSet
{
  private:

  const Context& ctxt_;
  std::vector<Constraint *> constraint_list;

  public:

  ConstraintSet(const Context& ctxt) : ctxt_(ctxt) {}
  bool define_constraint(AtomSet &atoms, int argc, char **argv);
  bool set_constraint(int argc, char **argv);
  bool delete_constraint(int argc, char **argv);
  void list_constraints(std::ostream &os);
  int size(void) const { return constraint_list.size(); }
  void enforce(AtomSet& atoms);
  void enforce_r(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &rp);
  void enforce_v(const std::vector<std::vector<double> > &r0,
                 std::vector<std::vector<double> > &v0);
  void compute_forces(const std::vector<std::vector<double> > &r0,
                      const std::vector<std::vector<double> > &f);
  void update_constraints(double dt);
  void setup(AtomSet& atoms);
};
#endif
