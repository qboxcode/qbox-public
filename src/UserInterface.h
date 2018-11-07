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
// UserInterface.h:
//
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include <iostream>
#include <string>
#include <cstring> // strncpy(), strtok()
#include <cstdlib> // free(), system()
#include <list>

class UserInterface;

class Cmd
{
  public:
  UserInterface *ui;
  virtual const char *name(void) const = 0;
  virtual const char *help_msg(void) const = 0;
  virtual int action(int argc, char **argv) = 0;
  virtual ~Cmd(void) {}
};

class Var
{
  public:
  UserInterface *ui;
  virtual const char *name ( void ) const = 0;
  virtual int set ( int argc, char **argv ) = 0;
  virtual std::string print ( void ) const = 0;
  virtual ~Var(void) {}
};

class UserInterface
{
  private:

  int readCmd(std::string& s, std::istream &is, bool echo);
  void execCmd(std::string s, std::string prompt);
  bool terminate_;
  bool onpe0_;

  public:

  std::list<Cmd*> cmdlist;
  std::list<Var*> varlist;

  void addCmd(Cmd *newcmd)
  {
    newcmd->ui = this;
    cmdlist.push_back( newcmd );
  };

  Cmd *findCmd(const char *cmdname)
  {
    std::list<Cmd*>::iterator cmd;
    for ( cmd = cmdlist.begin();
          (cmd != cmdlist.end() && (strcmp((*cmd)->name(),cmdname)));
          cmd++ );

    if ( cmd != cmdlist.end() )
    {
      return (*cmd);
    }
    else
    {
      return 0;
    }
  };

  void addVar(Var *newvar)
  {
    newvar->ui = this;
    varlist.push_back( newvar );
  };

  Var *findVar(const char *varname)
  {
    std::list<Var*>::iterator var;
    for ( var = varlist.begin();
          (var != varlist.end() && (strcmp((*var)->name(),varname)));
          var++ );

    if ( var != varlist.end() )
    {
      return (*var);
    }
    else
    {
      return 0;
    }
  };

  void processCmds(std::istream &cmdstream, const std::string prompt,
                   bool echo);
  void processCmdsServer(std::string inputfilename, std::string outputfilename,
                         const std::string prompt, bool echo);

  void terminate(void) { terminate_ = true; }

  bool onpe0(void) const { return onpe0_; }

  UserInterface(void);
  ~UserInterface(void);
};
#endif
