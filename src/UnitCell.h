////////////////////////////////////////////////////////////////////////////////
//
// UnitCell.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: UnitCell.h,v 1.3 2003-09-08 23:28:32 fgygi Exp $

#ifndef UNITCELL_H
#define UNITCELL_H

#include "D3vector.h"
using namespace std;

class UnitCell
{
  private:

  D3vector a_[3];
  D3vector b_[3];
  double volume_;
  D3vector an_[13];
  D3vector bn_[13];
  double an2h_[13];
  double bn2h_[13];
  
  public:

  const D3vector& a(int i) const { return a_[i]; }
  const D3vector& b(int i) const { return b_[i]; }
  
  UnitCell(void) { set(D3vector(0,0,0),D3vector(0,0,0),D3vector(0,0,0)); }  
  explicit UnitCell(const D3vector& a0, const D3vector& a1, const D3vector& a2)
  { set(a0,a1,a2); }
  
  void set(const D3vector& a0, const D3vector& a1, const D3vector& a2);
  double volume(void) const { return volume_; }
  
  bool in_ws(const D3vector& v) const;
  void fold_in_ws(D3vector& v) const;
  bool in_bz(const D3vector& k) const;
  void fold_in_bz(D3vector& k) const;
  
  bool encloses(const UnitCell& c) const;
  bool contains(D3vector v) const;
  
  void print(ostream& os) const;  
  bool operator==(const UnitCell& c) const;
};
ostream& operator << ( ostream& os, const UnitCell& cell );
#endif
