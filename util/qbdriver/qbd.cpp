////////////////////////////////////////////////////////////////////////////////
//
// qb_driver.cpp
//
////////////////////////////////////////////////////////////////////////////////
//
// use: qb_driver qb_input_file qb_output_file
//
// qb_driver sends commands to the server, via the file qb_input_file
// It checks for the presence of a link named "qb_input_file.lock"
// before writing additional commands
// compile with: g++ -o qb_driver qb_driver.cpp

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

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

int main(int argc, char** argv)
{
  char* qb_infilename = argv[1];
  char* qb_outfilename = argv[2];
  string lockfilename = string(qb_infilename) + ".lock";

  ofstream qb_infile;
  ifstream qb_outfile;

  // send command to execute init.i script
  wait_for_file(lockfilename);
  qb_infile.open(qb_infilename,ios_base::trunc);
  qb_infile << "init.i" << endl;
  cout << " sent init.i cmd" << endl;
  qb_infile.close();
  sync();
  remove(lockfilename.c_str());
  cout << " lock file removed" << endl;
  usleep(500000);

  // extract element <etotal> from output
  string element_name("etotal");
  double etotal;

  bool done = false;
  //while ( !done )
  for ( int iter = 0; iter < 5; iter++ )
  {
    cout << " loop start" << endl;

    wait_for_file(lockfilename);

    cout << " reading output ... ";
    qb_outfile.open(qb_outfilename);
    // get output from server
    // read output until EOF
    string s, qb_output;
    //qb_outfile.sync();
    qb_outfile.clear();
    while ( qb_outfile )
    {
      getline(qb_outfile,s);
      qb_output += s;
      qb_output += '\n';
    }
    cout << qb_output << endl;
    cout << "done" << endl;
    qb_outfile.close();

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
      cout << " element " << element_name << " not found in qb_output" << endl;
    }

    //   analyze data
    //
    // generate a random move in a cube of side 0.04

    const double amplitude = 0.02;
    double dx = amplitude * (2.0*drand48()-1.0);
    double dy = amplitude * (2.0*drand48()-1.0);
    double dz = amplitude * (2.0*drand48()-1.0);

    // prepare next commands
    // send next command to server

    qb_infile.open(qb_infilename,ios_base::trunc);
    qb_infile << "move C by " << dx << " " << dy << " " << dz << endl;
    qb_infile << "run 0 30" << endl;
    qb_infile.close();
    sync();
    remove(lockfilename.c_str());
    cout << " lock file  removed" << endl;
    usleep(100000);
  }

  wait_for_file(lockfilename);
  qb_infile.open(qb_infilename,ios_base::trunc);
  qb_infile << "quit" << endl;
  qb_infile.close();
  sync();
  cout << " sent quit cmd" << endl;
  remove(lockfilename.c_str());
  cout << " lock file removed" << endl;

  return 0;
}
