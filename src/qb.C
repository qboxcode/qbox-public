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
// qb.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: qb.C,v 1.60 2008-09-08 15:56:20 fgygi Exp $

#include <iostream>
#include <string>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>
#if AIX
#include<filehdr.h>
#endif
#ifdef USE_APC
#include "apc.h"
#endif

#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"

#include "Context.h"
#include "UserInterface.h"
#include "Sample.h"
#include "Timer.h"

#include "AngleCmd.h"
#include "AtomCmd.h"
#include "ComputeMLWFCmd.h"
#include "ConstraintCmd.h"
#include "DistanceCmd.h"
#include "FoldInWsCmd.h"
#include "HelpCmd.h"
#include "KpointCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "MoveCmd.h"
#include "PrintCmd.h"
#include "QuitCmd.h"
#include "RandomizeWfCmd.h"
#include "ResetVcmCmd.h"
#include "RunCmd.h"
#include "SaveCmd.h"
#include "SetCmd.h"
#include "SpeciesCmd.h"
#include "StatusCmd.h"
#include "TorsionCmd.h"

#include "AtomsDyn.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "ChargeMixCoeff.h"
#include "ChargeMixRcut.h"
#include "Debug.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecuts.h"
#include "Emass.h"
#include "ExtStress.h"
#include "FermiTemp.h"
#include "Dt.h"
#include "Nempty.h"
#include "NetCharge.h"
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

#if BGLDEBUG
#include <rts.h>
#endif

int main(int argc, char **argv, char **envp)
{
  Timer tm;
  tm.start();

#if USE_MPI
  MPI_Init(&argc,&argv);
#endif
#if USE_APC
  ApcInit();
#endif

#if BGLDEBUG
  {
    int myrank,mysize;
    BGLPersonality personality;
    rts_get_personality (&personality, sizeof(personality));

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    cout << myrank << ": at "
         << personality.xCoord << " "
         << personality.yCoord << " "
         << personality.zCoord << endl;
  }
#endif

  {
  Context ctxt;

  if ( ctxt.onpe0() )
  {
  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:simulation xmlns:fpmd=\"" << qbox_xmlns() << "\">" << endl;
  cout << "\n";
  cout << "                   ============================\n";
  cout << "                   I qbox "
       << setw(17) << left << release() << "   I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I                          I\n";
  cout << "                   I http://eslab.ucdavis.edu I\n";
  cout << "                   ============================\n\n";
  cout << "\n";
  cout << "<release> " << release() << " " << TARGET << " </release>" << endl;

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
  // remove angle brackets from processor name for XML compatibility
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
  {
    if ( processor_name[i] == '<' ) processor_name[i] = '(';
    if ( processor_name[i] == '>' ) processor_name[i] = ')';
  }
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
    {
      cout << "<process id=\"" << ip << "\"> " << buf
           << " </process>" << endl;
    }
  }
  if ( ctxt.onpe0() )
    cout << "</mpi_processes>" << endl;
#endif // USE_MPI

  Sample* s = new Sample(ctxt);

  UserInterface ui;

  ui.addCmd(new AngleCmd(s));
  ui.addCmd(new AtomCmd(s));
  ui.addCmd(new ComputeMLWFCmd(s));
  ui.addCmd(new ConstraintCmd(s));
  ui.addCmd(new DistanceCmd(s));
  ui.addCmd(new FoldInWsCmd(s));
  ui.addCmd(new HelpCmd(s));
  ui.addCmd(new KpointCmd(s));
  ui.addCmd(new ListAtomsCmd(s));
  ui.addCmd(new ListSpeciesCmd(s));
  ui.addCmd(new LoadCmd(s));
  ui.addCmd(new MoveCmd(s));
  ui.addCmd(new PrintCmd(s));
  ui.addCmd(new QuitCmd(s));
  ui.addCmd(new RandomizeWfCmd(s));
  ui.addCmd(new ResetVcmCmd(s));
  ui.addCmd(new RunCmd(s));
  ui.addCmd(new SaveCmd(s));
  ui.addCmd(new SetCmd(s));
  ui.addCmd(new SpeciesCmd(s));
  ui.addCmd(new StatusCmd(s));
  ui.addCmd(new TorsionCmd(s));

  ui.addVar(new AtomsDyn(s));
  ui.addVar(new Cell(s));
  ui.addVar(new CellDyn(s));
  ui.addVar(new CellLock(s));
  ui.addVar(new ChargeMixCoeff(s));
  ui.addVar(new ChargeMixRcut(s));
  ui.addVar(new CellMass(s));
  ui.addVar(new Debug(s));
  ui.addVar(new Ecut(s));
  ui.addVar(new Ecutprec(s));
  ui.addVar(new Ecuts(s));
  ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new FermiTemp(s));
  ui.addVar(new Dt(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new NetCharge(s));
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

  if ( argc == 2 )
  {
    // input file was given as a command line argument
    bool echo = true;
    ifstream in;
    if ( ctxt.onpe0() )
    {
      in.open(argv[1],ios::in);
    }
    ui.processCmds(in, "[qbox]", echo);
  }
  else
  {
    // use standard input
    bool echo = !isatty(0);
    ui.processCmds(cin, "[qbox]", echo);
  }

  // exit using the quit command when a encountering EOF in a script
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( ctxt.onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</fpmd:simulation>" << endl;
  }

  } // end of Context scope
#if USE_APC
  ApcFinalize();
#endif
#if USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
