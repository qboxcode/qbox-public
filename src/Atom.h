////////////////////////////////////////////////////////////////////////////////
//
// Atom.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Atom.h,v 1.2 2003-05-16 16:14:00 fgygi Exp $

#ifndef ATOM_H
#define ATOM_H

#include "D3vector.h"
#include <string>
using namespace std;

class Atom
{
  private:
  
  string name_;
  string species_;
  D3vector position_;
  D3vector velocity_;

  public:

  Atom (string name, string species, D3vector position, D3vector velocity);
  string name(void) { return name_; };
  string species(void) { return species_; };
  D3vector position(void) { return position_; };
  D3vector velocity(void) { return velocity_; };
  void set_position(D3vector p) { position_ = p; };
  void set_velocity(D3vector v) { velocity_ = v; };
  void block(void) { velocity_ = D3vector(0.0,0.0,0.0); };
};

ostream& operator << ( ostream &os, Atom &a );
#endif
