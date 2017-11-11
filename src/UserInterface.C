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
// UserInterface.C: definition of readCmd and processCmds
//
////////////////////////////////////////////////////////////////////////////////

#include "UserInterface.h"
#include "qbox_xmlns.h"
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <unistd.h> // fsync()
#include <stdio.h> // fopen(), fclose(), fprintf()
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
int UserInterface::readCmd(char *s, int max, istream &fp, bool echo)
{
  int ch, i = 0;
  while ( (ch = fp.get()) != EOF &&  !( ch == '\n' || ch ==';' || ch == '#') )
  {
    if ( ch == '\\' ) // line continuation character
    {
      // check if backslash is followed by a newline
      ch = fp.get();
      if ( ch == '\n' )
      {
        // backslash followed by newline, do nothing
      }
      else
      {
        // backslash not followed by newline
        if ( i < max - 1 )
          s[i++] = '\\';
        if ( i < max - 1 )
          s[i++] = ch;
      }
    }
    if (i < max - 1)
      s[i++] = ch;
  }

  if (max > 0) s[i] = '\0';  /* add terminating NULL */

  if ( fp.eof() )
    return 0;             /* return 0 for end of file */

  // output command line if reading from a script
  if ( echo && i > 0 ) cout << "<cmd>" << s << "</cmd>" << endl;

  if ( ch == '#' )
  {
    if ( echo ) cout << "<cmd>#";
    while ( (ch = fp.get()) != EOF && !( ch == '\n' ) )
    {
      if ( echo ) cout << (char) ch;
    }
    if ( echo && ch=='\n' ) cout << "</cmd>" << endl;
    if ( !(ch == '\n') )
      return 0;             /* return 0 for end of file */
  }

  return 1; // a command was read
}

////////////////////////////////////////////////////////////////////////////////
void UserInterface::processCmds ( istream &cmdstream, const char *prompt,
  bool echo)
{
#if DEBUG
  cout << "UserInterface::processCmds: prompt="
       << prompt << " echo=" << echo << endl;
#endif
  // read and process commands from cmdstream until done
  const int cmdlinemax = 1024;
  char cmdline[cmdlinemax];
  list<Cmd*>::iterator cmd;
  const char *separators = " ;\t";
  int done=0,cmd_read,status;

  if ( onpe0_ )
    cout << prompt << " ";

  while ( !done )
  {
    if ( onpe0_ )
    {
      for ( int i = 0; i < cmdlinemax; i++ )
        cmdline[i] = '\0';
      // readCmd returns 1 if a command is read, 0 if at EOF
      cmd_read = readCmd(cmdline, cmdlinemax, cmdstream, echo );
      done = !cmd_read;
    }
    MPI_Bcast(&cmdline[0],cmdlinemax,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&cmd_read,1,MPI_INT,0,MPI_COMM_WORLD);

    if ( cmd_read )
    {
      // cmdline contains a string of tokens terminated by '\0'
      // cout << " command line is: " << cmdline << endl;

      // comment lines: start with '#'
      int i;
      if ( cmdline[i=strspn(cmdline," ")] == '#' )
      {
        // cout << " comment line" << endl;
        // do nothing, not even write prompt
      }
      else if ( cmdline[i=strspn(cmdline," ")] == '!' )
      {
        // shell escape commands start with '!'
        // cout << " shell escape" << endl;
        if ( onpe0_ )
        {
          system ( &cmdline[i+1] );
          cout << prompt << " ";
        }
      }
      else
      {
        // cout << " command split in the following tokens:" << endl;

        // scan tokens and build argument list
        list<char*> arglist;
        int ntok = 0;
        char* tok = strtok(cmdline, separators);
        while ( tok != 0 )
        {
          arglist.push_back(tok);
          ntok++;
          // cout << "\"" << tok << "\"" << endl;
          tok = strtok(0,separators);
        }
        // cout << " total of " << ntok << " tokens" << endl;
        // arglist.dump();

        // build ac and av
        int ac = ntok;
        char **av = new char *[ntok+1];
        i = 0;
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
            if ( onpe0_ )
            {
              cmdstr.open(av[0],ios::in);
              status = !cmdstr;
            }
            MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
            if ( !status )
            {
              // create new prompt in the form: prompt<filename>
              char *newprompt=0;
              if ( onpe0_ )
              {
                newprompt = new char[strlen(prompt)+strlen(av[0])+4];
                newprompt = strcpy(newprompt,prompt);
                newprompt = strcat(newprompt,"[");
                newprompt = strcat(newprompt,av[0]);
                newprompt = strcat(newprompt,"]");
              }
                // MPI: process commands on all processes.
                // Note: newprompt is 0 on all processes > 0
                // Note: echo == true for scripts
              processCmds (cmdstr, newprompt, true);
              if ( onpe0_ )
                delete newprompt;
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
      }
    }
    else // if cmd_read
    {
      // found eof
      done = true;
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
  string outputfilename, const char *prompt, bool echo)
{
  // process commands in server mode
  // read commands from inputfilename
  // redirect cout to outputfilename

  char cmdline[256];
  list<Cmd*>::iterator cmd;
  const char *separators = " ;\t";
  int i,done=0,cmd_read,status;

  string lockfilename = inputfilename + ".lock";

  // in server mode, redirect output to stream qbout
  streambuf *cout_buf;
  ifstream qbin;

  ofstream tstfile;
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
        cmd_read = readCmd(cmdline, 256, qbin, echo );
        cout << prompt << " " << cmdline << endl;
      }
      MPI_Bcast(&cmdline[0],256,MPI_CHAR,0,MPI_COMM_WORLD);
      MPI_Bcast(&cmd_read,1,MPI_INT,0,MPI_COMM_WORLD);

      if ( cmd_read )
      {
        if ( onpe0_ )
        {
          cerr << " read cmd: " << cmdline << endl;
          cerr << " executing ... ";
        }

        // cmdline contains a string of tokens terminated by '\0'
        // cout << " command line is: " << cmdline << endl;

        // comment lines: start with '#'
        if ( cmdline[i=strspn(cmdline," ")] == '#' )
        {
          // cout << " comment line" << endl;
          // do nothing, not even write prompt
        }
        else if ( cmdline[i=strspn(cmdline," ")] == '!' )
        {
          // shell escape commands start with '!'
          // cout << " shell escape" << endl;
          if ( onpe0_ )
          {
            system ( &cmdline[i+1] );
            // cout << prompt << " ";
          }
        }
        else
        {
          // cout << " command split in the following tokens:" << endl;

          // scan tokens and build argument list
          list<char*> arglist;
          int ntok = 0;
          char* tok = strtok(cmdline, separators);
          while ( tok != 0 )
          {
            arglist.push_back(tok);
            ntok++;
            // cout << "\"" << tok << "\"" << endl;
            tok = strtok(0,separators);
          }
          // cout << " total of " << ntok << " tokens" << endl;
          // arglist.dump();

          // build ac and av
          int ac = ntok;
          char **av = new char *[ntok+1];
          i = 0;
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
              cerr << " execute command " << cmdptr->name() << endl;
#endif
              cmdptr->action(ac,av);
              MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
              cerr << " command completed " << cmdptr->name() << endl;
#endif
            }
            else
            {
              // command is not in the command list, check for script files
              ifstream cmdstr;
              if ( onpe0_ )
              {
                cmdstr.open(av[0],ios::in);
                status = !cmdstr;
              }
              MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
              if ( !status )
              {
                // create new prompt in the form: prompt<filename>
                char *newprompt=0;
                if ( onpe0_ )
                {
                  newprompt = new char[strlen(prompt)+strlen(av[0])+4];
                  newprompt = strcpy(newprompt,prompt);
                  newprompt = strcat(newprompt,"[");
                  newprompt = strcat(newprompt,av[0]);
                  newprompt = strcat(newprompt,"]");
                }
                // MPI: process commands on all processes.
                // Note: newprompt is 0 on all processes > 0
                // Note: echo == true for scripts
                processCmds (cmdstr, newprompt, true);
                if ( onpe0_ )
                  delete newprompt;
              }
              else
              {
                if ( onpe0_ )
                  cout << " No such command or file name: " << tok << endl;
              }
            }
          }
          delete [] av;
        }

        if ( onpe0_ )
          cerr << "done" << endl;

        // check if terminate_ flag was set during command execution
        if ( onpe0_ )
          done = terminate_;
        MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);

      } // if cmd_read

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
