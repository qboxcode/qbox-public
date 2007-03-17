////////////////////////////////////////////////////////////////////////////////
//
// NameOf.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NameOf.h,v 1.2 2007-03-17 01:14:00 fgygi Exp $

#ifndef NAMEOF_H
#define NAMEOF_H

#include <string>

// predicate class for searching T* containers by name
// T must be a pointer type to something that has a name() member
template <class T>
class NameOf
{
  public:
  std::string name;

  NameOf<T>(std::string s) : name(s) {};
  bool operator() (T t) const
  {
    return t->name() == name;
  }
};
#endif
