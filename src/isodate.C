////////////////////////////////////////////////////////////////////////////////
//
// isodate.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: isodate.C,v 1.1 2004-05-20 00:22:20 fgygi Exp $

#include "isodate.h"
#include <ctime>
std::string isodate(void)
{
  const time_t t = time(NULL);
  struct tm* tms = gmtime(&t);
  char s[32];
  const char* fmt = "%Y-%m-%dT%TZ";
  strftime(s,32,fmt,tms);
  return std::string(s);
}
