#include <uuid/uuid.h>
#include "uuid.h"
using namespace std;

std::string uuid(void)
{
  unsigned char uuid[16];
  char s[36];

  // generate uuid
  uuid_generate_time(uuid);

  // convert to string
  uuid_unparse(uuid,s);

  return string(s);
}
