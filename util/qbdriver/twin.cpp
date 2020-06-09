////////////////////////////////////////////////////////////////////////////////
//
// twin.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// use: twin qbin qbout
//
// twin sends commands to two Qbox servers, via the files qbin_0, qbin_1
// It checks for the presence of a link named "qbin_<n>.lock"
// before writing additional commands
// compile with: g++ -o twin twin.cpp

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h> // stat()
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

void sendCmd(const string filename, const string cmd);
void wait_for_file(const string& lockfilename);
void wait_for_nofile(const string& lockfilename);

int main(int argc, char** argv)
{
  char* qb_infilename_root = argv[1];
  char* qb_outfilename_root = argv[2];

  // number of servers
  const int ns = 2;
  string qb_infilename[ns];
  string qb_outfilename[ns];
  string lockfilename[ns];
  for ( int i = 0; i < ns; i++ )
  {
    ostringstream i_str;
    i_str << i;
    qb_infilename[i] = string(qb_infilename_root) + "_" + i_str.str();
    qb_outfilename[i] = string(qb_outfilename_root) + "_" + i_str.str();
    lockfilename[i] = string(qb_infilename[i]) + ".lock";

    cout << " qb input file:  " << qb_infilename[i] << endl;
    cout << " qb output file: " << qb_outfilename[i] << endl;
    cout << " qb lock file:   " << lockfilename[i] << endl;
  }

  ifstream qb_outfile[2];

  // send commands to servers to execute the init.i script
  for ( int i = 0; i < ns; i++ )
  {
    wait_for_file(lockfilename[i]);
    sendCmd(qb_infilename[i],string("init.i\n"));
    remove(lockfilename[i].c_str());
  }
  usleep(100000);

  // extract element <etotal> from output
  string element_name("etotal");
  double etotal;

  for ( int iter = 0; iter < 5; iter++ )
  {
    // loop over servers
    for ( int i = 0; i < ns; i++ )
    {
      wait_for_file(lockfilename[i]);

      cout << " reading output from server " << i << " ... " << endl;
      qb_outfile[i].open(qb_outfilename[i].c_str());
      // get output from server
      // read output until EOF
      string s, qb_output;
      //qb_outfile.sync();
      qb_outfile[i].clear();
      while ( qb_outfile[i] )
      {
        getline(qb_outfile[i],s);
        qb_output += s;
        qb_output += '\n';
      }
      cout << i << ": " << qb_output << endl;
      cout << i << ": " << "done" << endl;
      qb_outfile[i].close();

      // parse output
      string start_tag = "<" + element_name + ">";
      string end_tag = "</" + element_name + ">";
      size_t pstart = qb_output.find(start_tag);
      size_t pend = qb_output.find(end_tag);
      size_t len;
      if ( pend != string::npos )
      {
        string buf;
        // element was found
        len = pend + end_tag.size() - pstart;
        buf = qb_output.substr(pstart,len);
        cout << " buf= " << buf << endl;
        istringstream is(buf);
        string dum;
        is >> dum >> etotal;
        cout << " etotal = " << etotal << endl;
      }
      else
      {
        cout << " element " << element_name
             << " not found in qb_output[" << i << "]" << endl;
      }

      // analyze data
      //
      // generate a random move in a cube of side 0.04

      const double amplitude = 0.02;
      double dx = amplitude * (2.0*drand48()-1.0);
      double dy = amplitude * (2.0*drand48()-1.0);
      double dz = amplitude * (2.0*drand48()-1.0);

      // prepare next commands
      // send next commands to server

      // write all commands into ostringstream os
      ostringstream os;
      os << "move C by " << dx << " " << dy << " " << dz << endl;
      os << "run 0 30" << endl;

      // write ostringstream to file
      sendCmd(qb_infilename[i],os.str());
      remove(lockfilename[i].c_str());
    }
    usleep(100000);
  }

  // send quit command to all servers
  for ( int i = 0; i < ns; i++ )
  {
    wait_for_file(lockfilename[i]);
    sendCmd(qb_infilename[i],string("quit\n"));
    remove(lockfilename[i].c_str());
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
void sendCmd(const string filename, const string str)
{
  FILE *fp = fopen(filename.c_str(),"w");
  fprintf(fp,"%s",str.c_str());
  fclose(fp);
  fsync(fileno(fp));
}

////////////////////////////////////////////////////////////////////////////////
void wait_for_file(const string& lockfilename)
{
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat(lockfilename.c_str(),&statbuf);
    usleep(100000);
  }
  while ( status != 0 );
}

////////////////////////////////////////////////////////////////////////////////
void wait_for_nofile(const string& lockfilename)
{
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat(lockfilename.c_str(),&statbuf);
    usleep(100000);
  }
  while ( status == 0 );
}
