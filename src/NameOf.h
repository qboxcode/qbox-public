////////////////////////////////////////////////////////////////////////////////
//
// NameOf.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: NameOf.h,v 1.1.1.1 2002-09-27 00:08:39 fgygi Exp $

#ifndef NAMEOF_H
#define NAMEOF_H

#include <string>
using namespace std;

// predicate class for searching T* containers by name
// T must be a pointer type to something that has a name() member
template <class T>
class NameOf
{
  public:
  string name;

  NameOf<T>(string s) : name(s) {};
  bool operator() (T t) const
  {
    return t->name() == name;
  }
};
#endif
