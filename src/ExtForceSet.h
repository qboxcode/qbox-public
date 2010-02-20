////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// ExtForceSet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExtForceSet.h,v 1.1 2010-02-20 23:13:02 fgygi Exp $

#ifndef EXTFORCESET_H
#define EXTFORCESET_H

#include <vector>
#include <string>

class Atom;
class AtomSet;
class ExtForce;
class Context;

class ExtForceSet
{
  private:

  const Context& ctxt_;
  std::vector<ExtForce *> extforce_list;

  public:

  ExtForceSet(const Context& ctxt) : ctxt_(ctxt) {}
  ~ExtForceSet();
  bool define_extforce(AtomSet &atoms, int argc, char **argv);
  bool set_extforce(int argc, char **argv);
  bool delete_extforce(int argc, char **argv);
  void list_extforces(std::ostream &os);
  int size(void) const { return extforce_list.size(); }
  double energy(const std::vector<std::vector<double> > &r,
                std::vector<std::vector<double> > &f) const;
  void setup(AtomSet& atoms);
  void reset(void);
};
#endif
