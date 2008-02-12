////////////////////////////////////////////////////////////////////////////////
//
// UserInterface.C: definition of readCmd and processCmds
//
////////////////////////////////////////////////////////////////////////////////
// $Id: UserInterface.C,v 1.6 2008-02-12 05:39:18 fgygi Exp $

#include "UserInterface.h"
#include <string>
#include <list>
#include <unistd.h> // isatty
#include <fstream>

#if USE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif
using namespace std;

////////////////////////////////////////////////////////////////////////////////
UserInterface::UserInterface(void) : terminate_(false)
{
#if USE_MPI
  int mype;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  onpe0_ = ( mype == 0 );
#else
  onpe0_ = true;
#endif
}

////////////////////////////////////////////////////////////////////////////////
char *UserInterface::readCmd(char *s, int max, istream &fp, bool echo)
{
  int ch, i = 0;
  while ( (ch = fp.get()) != EOF && !( ch == '\n' || ch ==';' || ch == '#') )
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
    else
    {
      if (i < max - 1)
        s[i++] = ch;
    }
  }
  if (max > 0) s[i] = '\0';  /* add terminating NULL */

  if ( !(ch == '\n' || ch == ';' || ch == '#') )
    return NULL;             /* return NULL for end of file */

  // output command line if reading from a script
  if ( echo ) cout << s;

  if ( ch == '#' )
  {
    if ( echo ) cout << '#';
    while ( (ch = fp.get()) != EOF && !( ch == '\n' ) )
    {
      if ( echo ) cout << (char) ch;
    }
    if ( !(ch == '\n') )
      return NULL;             /* return NULL for end of file */
  }

  return s;
}

////////////////////////////////////////////////////////////////////////////////
void UserInterface::processCmds ( istream &cmdstream, char *prompt, bool echo )
{
  // read and process commands from cmdstream until end of file is reached

  char cmdline[256];
  list<Cmd*>::iterator cmd;
  char *tok;
  const char *separators = " ;\t";
  int i,done,status;

  if ( onpe0_ )
    cout << prompt << " ";

  // read a command terminated by '\n' or ';'
  if ( onpe0_ )
  {
    done = terminate_ || !readCmd(cmdline, 256, cmdstream, echo );
    // readCmd returns cmdline if a command is read, NULL at EOF
  }
#if USE_MPI
  MPI_Bcast(&cmdline[0],256,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  while ( !done )
  {
    if ( onpe0_ )
      cout << endl;

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
        cout << prompt << " ";
      }
    }
    else
    {
      // cout << " command split in the following tokens:" << endl;

      // scan tokens and build argument list
      list<char*> arglist;
      int ntok = 0;
      tok = strtok(cmdline, separators);
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
#if USE_MPI
          MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
          cout << " execute command " << cmdptr->name() << endl;
#endif
#endif
          cmdptr->action(ac,av);
#if USE_MPI
          MPI_Barrier(MPI_COMM_WORLD);
#if DEBUG
          cout << " command completed" << cmdptr->name() << endl;
#endif
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
#if USE_MPI
          MPI_Bcast(&status,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
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

    // read a command terminated by '\n' or ';'
    if ( onpe0_ )
    {
      done = terminate_ || !readCmd(cmdline, 256, cmdstream, echo );
    }
#if USE_MPI
    MPI_Bcast(&cmdline[0],256,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

  }
  if ( onpe0_ )
    cout << " End of command stream " << endl;
}
