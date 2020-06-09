//
// testreassign.cpp
//
// test the functionality of:
// - reassign the streambuf of std::cout to an ostringstream
// - write and sync using C library functions
//
// This program tests the functionality needed in UserInterface when
// operating in server mode. In that case, the std::cout stream is redirected
// to an ostringstream. The contents of the ostringstream are written at the
// end using C library functions, which allow for the use of the fsync() call.
//

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
int main()
{
  // write something to std::cout before reassign
  std::cout << "initial write on std::cout" << std::endl;

  // streambuf pointers
  std::streambuf *qbout_buf;
  std::streambuf *cout_buf;

  // save copy of cout streambuf
  cout_buf = std::cout.rdbuf();
  // create an ostringstream
  std::ostringstream os;
  qbout_buf = os.rdbuf();
  // redirect std::cout to os
  std::cout.rdbuf(qbout_buf);

  // the following output should go to the ostringstream
  std::cout << " output written to cout after redirect" << std::endl;
  std::cout.flush();

  // write contents of os to file "out.txt"
  FILE *qboutfile = fopen("out.txt","w");
  fprintf(qboutfile,"%s",os.str().c_str());
  fclose(qboutfile);
  fsync(fileno(qboutfile));

  // restore cout streambuf
  std::cout.rdbuf(cout_buf);

  // write more output to std::cout
  std::cout << "more output on std::cout" << std::endl;
}
