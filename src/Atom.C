////////////////////////////////////////////////////////////////////////////////
//
// Atom.C:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Atom.C,v 1.3 2003-05-16 16:14:00 fgygi Exp $

#include "Atom.h"
#include <iomanip>
using namespace std;

Atom::Atom (string newname, string newspecies, D3vector pos, D3vector vel)
{
  name_ = newname;
  species_ = newspecies;
  position_ = pos;
  velocity_ = vel;
}

ostream& operator << ( ostream &os, Atom &a )
{
  os.setf(ios::left,ios::adjustfield);
  os << "  <atom name=\"" << a.name() << "\""
     << " species=\"" << a.species() << "\">\n"
     << "    <position> ";
  os.setf(ios::fixed,ios::floatfield);
  os.setf(ios::right,ios::adjustfield);
  os << setw(12) << setprecision(8) << a.position().x << " "
     << setw(12) << setprecision(8) << a.position().y << " "
     << setw(12) << setprecision(8) << a.position().z << "  "
     << " </position>\n"
     << "    <velocity> ";
  os.setf(ios::scientific,ios::floatfield);
  os << setw(13) << setprecision(6) << a.velocity().x << " "
     << setw(13) << setprecision(6) << a.velocity().y << " "
     << setw(13) << setprecision(6) << a.velocity().z
     << " </velocity>\n  </atom>" << endl;
  return os;
}
