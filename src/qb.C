////////////////////////////////////////////////////////////////////////////////
//
// qb.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: qb.C,v 1.33 2004-04-20 22:18:05 fgygi Exp $

const char* const release = "1.14.1a";
const char* const xmlns_url = "http://www.llnl.gov/casc/fpmd/qbox/1.0";

#include <iostream>
#include <string>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <ctime>
#include <cstdlib>
#if AIX || OSF1
#include<filehdr.h>
#endif

#include "Context.h"
#include "UserInterface.h"
#include "Sample.h"
#include "Timer.h"

#include "AtomCmd.h"
#include "HelpCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "PrintCmd.h"
#include "QuitCmd.h"
#include "RandomizeWfCmd.h"
#include "RunCmd.h"
#include "SaveCmd.h"
#include "SetCmd.h"
#include "SpeciesCmd.h"
#include "StatusCmd.h"

#include "AtomsDyn.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "Debug.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecuts.h"
#include "Emass.h"
#include "ExtStress.h"
#include "Dt.h"
#include "Nempty.h"
#include "Nrowmax.h"
#include "RefCell.h"
#include "Stress.h"
#include "Thermostat.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "Xc.h"

string isodate(void)
{
  const time_t t = time(NULL);
  struct tm* tms = gmtime(&t);
  char s[32];
  const char* fmt = "%Y-%m-%dT%TZ";
  strftime(s,32,fmt,tms);
  string st(s);
  return st;
}

int main(int argc, char **argv, char **envp)
{
  Timer tm;
  tm.start();

#if USE_MPI
  MPI_Init(&argc,&argv);
#endif

  {
  Context ctxt;
  
  if ( ctxt.onpe0() )
  {
  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<!--\n\n";
  cout << "                   ===========================\n";
  cout << "                   I qbox " 
       << setw(17) << left << release << "  I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I                         I\n";
  cout << "                   I            F.Gygi, LLNL I\n";
  cout << "                   I Copyright (c) 2003-2004 I\n";     
  cout << "                   ===========================\n\n";
  cout << "-->\n";
  cout << "<qbox:simulation xmlns:qbox=\"" << xmlns_url << "\">" << endl;
  cout << "<release> " << release << " " << TARGET << " </release>" << endl;

  // Identify executable name, checksum, size and link date
  if ( getlogin() != 0 ) 
    cout << "<user> " << getlogin() << " </user>" << endl;
#if AIX || OSF1
  // read filehdr for link time
  filehdr hdr;
  FILE *fx = fopen(argv[0],"r");
  if ( fx != 0 )
  {
    size_t sz = fread((void*)&hdr,sizeof(filehdr),1,fx);
    fclose(fx);
    string s = ctime((time_t*)&hdr.f_timdat);
    cout << "<linktime> " << s << " </linktime>" << endl;
  }
#endif
  
  // Identify platform
  {
    struct utsname un;
    uname (&un);
    cout << "<sysname> " << un.sysname << " </sysname>" << endl;
    cout << "<nodename> " << un.nodename << " </nodename>" << endl;
  }

  cout << "<start_time> " << isodate() << " </start_time>" << endl;

  }

#if USE_MPI
  // Print list of node names  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  char buf[MPI_MAX_PROCESSOR_NAME];
  int namelen;
  PMPI_Get_processor_name(processor_name,&namelen);
  if ( ctxt.onpe0() )
  {
    cout << "<mpi_processes count=\"" << ctxt.size() << "\">" << endl;
    cout << "<process id=\"" << ctxt.mype() << "\"> " << processor_name 
         << " </process>" << endl;
  }
  for ( int ip = 1; ip < ctxt.size(); ip++ )
  {
    MPI_Barrier(ctxt.comm());
    if ( ctxt.onpe0() )
    {
      MPI_Status status;
      MPI_Recv(&buf[0],MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
                   ip,ip,ctxt.comm(),&status);
    }
    else if ( ip == ctxt.mype() )
    {
      // send processor name to pe0
      MPI_Send(&processor_name[0],MPI_MAX_PROCESSOR_NAME,
        MPI_CHAR,0,ctxt.mype(),ctxt.comm());
    }
  if ( ctxt.onpe0() )
    cout << "<process id=\"" << ip << "\"> " << buf 
         << " </process>" << endl;
  }
  if ( ctxt.onpe0() )
    cout << "</mpi_processes>" << endl;
#endif // USE_MPI

  Sample* s = new Sample(ctxt);
  
  UserInterface ui;
  
  ui.addCmd(new AtomCmd(s));
  ui.addCmd(new HelpCmd(s));
  ui.addCmd(new ListAtomsCmd(s));
  ui.addCmd(new ListSpeciesCmd(s));
  ui.addCmd(new LoadCmd(s));
  ui.addCmd(new PrintCmd(s));
  ui.addCmd(new QuitCmd(s));
  ui.addCmd(new RandomizeWfCmd(s));
  ui.addCmd(new RunCmd(s));
  ui.addCmd(new SaveCmd(s));
  ui.addCmd(new SetCmd(s));
  ui.addCmd(new SpeciesCmd(s));
  ui.addCmd(new StatusCmd(s));
  
  ui.addVar(new AtomsDyn(s));
  ui.addVar(new Cell(s));
  ui.addVar(new CellDyn(s));
  ui.addVar(new CellLock(s));
  ui.addVar(new CellMass(s));
  ui.addVar(new Debug(s));
  ui.addVar(new Ecut(s));
  ui.addVar(new Ecutprec(s));
  ui.addVar(new Ecuts(s));
  ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new Dt(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new Nrowmax(s));
  ui.addVar(new RefCell(s));
  ui.addVar(new Stress(s));
  ui.addVar(new Thermostat(s));
  ui.addVar(new ThTemp(s));
  ui.addVar(new ThTime(s));
  ui.addVar(new ThWidth(s));
  ui.addVar(new WfDiag(s));
  ui.addVar(new WfDyn(s));
  ui.addVar(new Xc(s));
  
  bool echo = !isatty(0);
  ui.processCmds(cin, "[qbox]", echo);
  // exit using the quit command when a encountering EOF in a script
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( ctxt.onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</qbox:simulation>" << endl;
  }

  } // end of Context scope
#if USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
