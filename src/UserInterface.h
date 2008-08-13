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
// $ Id: $

#ifndef USER_INTERFACE_H
#define USER_INTERFACE_H

#include <iostream>
#include <string>
#include <iomanip>
#include <list>
#include <algorithm>

class UserInterface;

class Cmd
{
  public:
  UserInterface *ui;
  virtual char *name(void) const = 0;
  virtual char *help_msg(void) const = 0;
  virtual int action(int argc, char **argv) = 0;
};

class Var
{
  public:
  UserInterface *ui;
  virtual char *name ( void ) const = 0;
  virtual int set ( int argc, char **argv ) = 0;
  virtual std::string print ( void ) const = 0;
};

class UserInterface
{
  private:

  char *readCmd(char *s, int max, std::istream &fp, bool echo);
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

  Cmd *findCmd(char *cmdname)
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

  Var *findVar(char *varname)
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

  void processCmds(std::istream &cmdstream, char *prompt, bool echo);

  void terminate(void) { terminate_ = true; }

  bool onpe0(void) const { return onpe0_; }

  UserInterface(void);
};
#endif
