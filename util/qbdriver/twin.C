////////////////////////////////////////////////////////////////////////////////
//
// qb_driver_twin.C
//
////////////////////////////////////////////////////////////////////////////////
//
// use: qb_driver qb_input qb_output
//
// qb_driver sends commands to the server, via the file qb_input_1, qb_input_2
// It checks for the presence of a link named "qb_input_<n>.lock"
// before writing additional commands
// compile with: g++ -o qb_driver qb_driver.C

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

void wait_for_file(string& lockfilename);
void wait_for_nofile(string& lockfilename);

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

  ofstream qb_infile[2];
  ifstream qb_outfile[2];

  // send commands to servers to execute the init.i script
  for ( int i = 0; i < ns; i++ )
  {
    wait_for_file(lockfilename[i]);
    qb_infile[i].open(qb_infilename[i].c_str(),ios_base::trunc);
    qb_infile[i] << "init.i" << endl;
    cout << " sent init.i cmd to server " << i << endl;
    qb_infile[i].close();
    sync();
    remove(lockfilename[i].c_str());
    cout << " lock file of server " << i << " removed" << endl;
  }
  usleep(500000);

  // extract element <etotal> from output
  string element_name("etotal");
  double etotal;

  for ( int iter = 0; iter < 5; iter++ )
  {
    cout << " loop start" << endl;

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
      // send next command to server

      qb_infile[i].open(qb_infilename[i].c_str(),ios_base::trunc);
      qb_infile[i] << "move C by " << dx << " " << dy << " " << dz << endl;
      qb_infile[i] << "run 0 30" << endl;
      qb_infile[i].close();
      sync();
      remove(lockfilename[i].c_str());
      cout << " lock file " << i << " removed" << endl;
    }
    usleep(100000);
  }

  for ( int i = 0; i < ns; i++ )
  {
    wait_for_file(lockfilename[i]);
    qb_infile[i].open(qb_infilename[i].c_str(),ios_base::trunc);
    qb_infile[i] << "quit" << endl;
    qb_infile[i].close();
    sync();
    cout << " sent quit cmd to server " << i << endl;
    remove(lockfilename[i].c_str());
    cout << " lock file " << i << " removed" << endl;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
void wait_for_file(string& lockfilename)
{
  cerr << " waiting for " << lockfilename << endl;
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat(lockfilename.c_str(),&statbuf);
    // cout << " status = " << status << endl;
    usleep(100000);
  }
  while ( status != 0 );
}

////////////////////////////////////////////////////////////////////////////////
void wait_for_nofile(string& lockfilename)
{
  cerr << " waiting for no " << lockfilename << endl;
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat(lockfilename.c_str(),&statbuf);
    // cout << " status = " << status << endl;
    usleep(100000);
  }
  while ( status == 0 );
}
