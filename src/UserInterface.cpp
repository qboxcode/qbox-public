////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// UserInterface.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "UserInterface.h"
#include "qbox_xmlns.h"
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <unistd.h> // fsync()
#include <cstdio> // fopen(), fclose(), fprintf()
#include <sys/types.h>
#include <sys/stat.h> // stat()
#include <mpi.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void wait_for_no_file(const string& lockfilename)
{
  //cerr << " waiting for no " << lockfilename << endl;
  struct stat statbuf;
  int status;
  do
  {
    // stat returns 0 if the file exists
    status = stat(lockfilename.c_str(),&statbuf);
    usleep(100000);
    // cerr << ".";
  }
  while ( status == 0 );
}

////////////////////////////////////////////////////////////////////////////////
UserInterface::UserInterface(void) : terminate_(false)
{
  int mype = 0;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  onpe0_ = ( mype == 0 );
}

////////////////////////////////////////////////////////////////////////////////
UserInterface::~UserInterface(void)
{
  std::list<Cmd*>::iterator cmd;
  for ( cmd = cmdlist.begin(); cmd != cmdlist.end(); cmd++ )
    delete *cmd;
  std::list<Var*>::iterator var;
  for ( var = varlist.begin(); var != varlist.end(); var++ )
    delete *var;
}

////////////////////////////////////////////////////////////////////////////////
int UserInterface::readCmd(string& s, istream &is, bool echo)
{
  getline(is,s);

  if ( is.eof() )
    return 0;

  // process line continuation if last character is '\'
  while ( (s.size() > 0) && (s[s.size()-1] == '\\') )
  {
    s[s.size()-1] = ' ';
    string stmp;
    getline(is,stmp);
    s += stmp;
  }

  if ( echo )
    cout << "<cmd>" << s << "</cmd>" << endl;

  // remove comments before processing
  string::size_type pos = s.find_first_of("#");
  if ( pos != string::npos )
    s = s.substr(0,pos);
#ifdef DEBUG
  cout << "readCmd: after removing comments: @" << s << "@" << endl;
#endif

  // remove leading white space
  pos = s.find_first_not_of(" \t");
  if ( pos != string::npos )
    s = s.substr(pos);
  else
    s.clear();
#ifdef DEBUG
  cout << "readCmd: after removing leading space: @" << s << "@" << endl;
#endif

  // remove trailing white space
  pos = s.find_last_not_of(" \t");
  if ( pos != string::npos )
  {
    pos++;
    s = s.substr(0,pos);
  }
#ifdef DEBUG
  cout << "readCmd: after removing trailing space: @" << s << "@" << endl;
#endif

  return 1; // a command was read into s
}

////////////////////////////////////////////////////////////////////////////////
void UserInterface::processCmds(istream &cmdstream, const string prompt,
  bool echo)
{
#if DEBUG
  cout << "UserInterface::processCmds: prompt="
       << prompt << " echo=" << echo << endl;
#endif
  // read and process commands from cmdstream until done
  string cmdline;
  list<Cmd*>::iterator cmd;
  int done = 0, cmd_read = 0;

  if ( onpe0_ )
    cout << prompt << " ";

  while ( !done )
  {
    if ( onpe0_ )
    {
      // readCmd returns 1 if a command is read, 0 if at EOF
      cmd_read = readCmd(cmdline, cmdstream, echo );
      done = !cmd_read;
    }

    MPI_Bcast(&cmd_read,1,MPI_INT,0,MPI_COMM_WORLD);

    if ( cmd_read )
    {
      int n = cmdline.size();
      MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
      char* buf = new char[n];
      if ( onpe0_ )
        strncpy(buf,cmdline.c_str(),n);
      MPI_Bcast(buf,n,MPI_CHAR,0,MPI_COMM_WORLD);
      cmdline = string(buf,n);
      delete [] buf;

      execCmd(cmdline,prompt);
    }

    if ( onpe0_ )
      done |= terminate_;

    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
  }

  if ( onpe0_ )
    cout << " End of command stream " << endl;
}

////////////////////////////////////////////////////////////////////////////////
void UserInterface::processCmdsServer ( string inputfilename,
  string outputfilename, const string prompt, bool echo)
{
  // process commands in server mode
  // read commands from inputfilename
  // redirect cout to outputfilename

  string cmdline;
  list<Cmd*>::iterator cmd;
  int done = 0, cmd_read = 0;

  string lockfilename = inputfilename + ".lock";

  // in server mode, redirect output to stream qbout
  streambuf *cout_buf;
  ifstream qbin;
  ostringstream os;

  while ( !done )
  {
    if ( onpe0_ )
    {
      // create file to signal that Qbox is waiting for a command on qbin
      FILE *lockfile = fopen(lockfilename.c_str(),"w");
      fprintf(lockfile,"1");
      fclose(lockfile);
      fsync(fileno(lockfile));
      usleep(100000);
      wait_for_no_file(lockfilename.c_str());

      qbin.open(inputfilename.c_str());
      qbin.sync();
      qbin.clear();

      // clear os
      os.str("");

      // save copy of cout streambuf
      cout_buf = cout.rdbuf();
      // redirect cout to os
      cout.rdbuf(os.rdbuf());
      // cerr << " processCmdsServer: cout streambuf redirected" << endl;

      cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
      cout << "<fpmd:simulation xmlns:fpmd=\"" << qbox_xmlns() << "\">" << endl;
    }

    do
    {
      if ( onpe0_ )
      {
        // readCmd returns 1 if a command is read, 0 if at EOF
        cmd_read = readCmd(cmdline, qbin, echo);
        cout << prompt << " " << cmdline << endl;
      }

      MPI_Bcast(&cmd_read,1,MPI_INT,0,MPI_COMM_WORLD);

      if ( cmd_read )
      {
        int n = cmdline.size();
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        char* buf = new char[n];
        if ( onpe0_ )
          strncpy(buf,cmdline.c_str(),n);
        MPI_Bcast(buf,n,MPI_CHAR,0,MPI_COMM_WORLD);
        cmdline = string(buf,n);
        delete [] buf;

        execCmd(cmdline,prompt);
      }

      if ( onpe0_ )
        done |= terminate_;

      MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);

    } while ( !done && cmd_read );

    if ( onpe0_ )
    {
      qbin.close();
      cout << " End of command stream " << endl;
      cout << "</fpmd:simulation>" << endl;
      cout.flush();

      // write ostringstream contents to output file
      FILE *qboutfile = fopen(outputfilename.c_str(),"w");
      fprintf(qboutfile,"%s",os.str().c_str());
      fclose(qboutfile);
      fsync(fileno(qboutfile));

      // restore cout streambuf
      cout.rdbuf(cout_buf);
      // cerr << " processCmdsServer: cout streambuf reassigned" << endl;
    }

    // wait before retrying
    usleep(200000);

  } // while !done

  if ( onpe0_ )
  {
    // remove lock file
    remove(lockfilename.c_str());
  }
}

////////////////////////////////////////////////////////////////////////////////
void UserInterface::execCmd(string cmdline, string prompt)
{
  // check if line contains only white space or comments
  if ( cmdline.size() == 0 )
  {
#ifdef DEBUG
    cout << " empty command line" << endl;
#endif
    if ( onpe0_ )
      cout << prompt << " ";
  }
  else if ( cmdline[0] == '!' )
  {
    // shell escape commands start with '!'
    if ( onpe0_ )
    {
      // remove '!' character
      cmdline[0] = ' ';
#ifdef DEBUG
      cout << " shell escape: " << cmdline.c_str() << endl;
#endif
      system( cmdline.c_str() );
      cout << prompt << " ";
    }
  }
  else
  {
#ifdef DEBUG
    cout << " command is: " << cmdline.c_str() << endl;
#endif
    // cout << " command split in the following tokens:" << endl;

    // scan tokens and build argument list
    list<char*> arglist;
    int ntok = 0;
    char* buf = strdup(cmdline.c_str());
    char* tok = strtok(buf, " \t");
    while ( tok != 0 )
    {
      arglist.push_back(tok);
      ntok++;
      // cout << "\"" << tok << "\"" << endl;
      tok = strtok(0," \t");
    }
    // cout << " total of " << ntok << " tokens" << endl;
    // arglist.dump();

    // build ac and av
    int ac = ntok;
    char **av = new char *[ntok+1];
    int i = 0;
    list<char*>::iterator iarg = arglist.begin();
    while ( iarg != arglist.end() )
    {
      av[i++] = *iarg++;
    }
    av[ntok] = 0;

    // write arguments
    // for ( i = 0; i < ac; i++ )
    // {
    //   cout << av[i] << endl;
    // }

    // search cmdlist for command

    tok = av[0];

    // check for empty command line
    if ( tok != 0 )
    {
      Cmd *cmdptr = findCmd(tok);

      if ( cmdptr )
      {
        MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
        cout << " execute command " << cmdptr->name() << endl;
#endif
        cmdptr->action(ac,av);
        MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
        cout << " command completed " << cmdptr->name() << endl;
#endif
      }
      else
      {
        // command is not in the command list, check for script files
        ifstream cmdstr;
        int status;
        if ( onpe0_ )
        {
          cmdstr.open(av[0],ios::in);
          status = !cmdstr;
        }
        MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
        if ( !status )
        {
          // create new prompt in the form: prompt<filename>
          string newprompt;
          if ( onpe0_ )
          {
            newprompt = prompt + "[" + av[0] + "]";
          }
            // MPI: process commands on all processes.
            // Note: newprompt is 0 on all processes > 0
            // Note: echo == true for scripts
          processCmds(cmdstr, newprompt, true);
        }
        else
        {
          if ( onpe0_ )
            cout << " No such command or file name: " << tok << endl;
        }
      }
    }
    delete [] av;
    if ( onpe0_ )
      cout << prompt << " ";
    free(buf);
  }
}
