//
// testlock.cpp
//
// test the creation and removal of a lock file using the C library
//

#include <stdio.h>    // fopen, fclose, fprintf
#include <sys/stat.h> // stat()
#include <unistd.h>   // fsync()
int main()
{
  // create a lock file "lock.dat"
  FILE *lockfile = fopen("lock.dat","w");
  fprintf(lockfile,"1");
  fclose(lockfile);
  fsync(fileno(lockfile));

  // test for the presence of the lock file
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat("lock.dat",&statbuf);
    usleep(100000);
  }
  while ( status != 0 );

  // remove the lock file
  remove("lock.dat");
  usleep(100000);

  // test for the absence of the file
  do
  {
    // stat returns 0 if the file exists
    status = stat("lock.dat",&statbuf);
    usleep(100000);
  }
  while ( status == 0 );
}
