////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2015 The Regents of the University of California
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

#include <iostream>
#include <string>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <cstdlib>
#include <fstream>

#include "isodate.h"
#include "release.h"
#include "qbox_xmlns.h"
#if USE_UUID
#include "uuid_str.h"
#endif
#ifdef _OPENMP
#include "omp.h"
#endif

#include "Context.h"
#include "UserInterface.h"
#include "Sample.h"
#include "Timer.h"

#include "AngleCmd.h"
#include "AtomCmd.h"
#include "ComputeMLWFCmd.h"
#include "ConstraintCmd.h"
#include "DistanceCmd.h"
#include "ExtForceCmd.h"
#include "FoldInWsCmd.h"
#include "HelpCmd.h"
#include "KpointCmd.h"
#include "ListAtomsCmd.h"
#include "ListSpeciesCmd.h"
#include "LoadCmd.h"
#include "MoveCmd.h"
#include "PartialChargeCmd.h"
#include "PlotCmd.h"
#include "PrintCmd.h"
#include "QuitCmd.h"
#include "RandomizeRCmd.h"
#include "RandomizeVCmd.h"
#include "RandomizeWfCmd.h"
#include "ResetRotationCmd.h"
#include "ResetVcmCmd.h"
#include "RescaleVCmd.h"
#include "ResponseCmd.h"
#include "RseedCmd.h"
#include "RunCmd.h"
#include "SaveCmd.h"
#include "SetCmd.h"
#include "SetVelocityCmd.h"
#include "SpeciesCmd.h"
#include "StatusCmd.h"
#include "StrainCmd.h"
#include "TorsionCmd.h"
#include "BisectionCmd.h"

#include "AlphaPBE0.h"
#include "AlphaRSH.h"
#include "AtomsDyn.h"
#include "BetaRSH.h"
#include "BlHF.h"
#include "BtHF.h"
#include "Cell.h"
#include "CellDyn.h"
#include "CellLock.h"
#include "CellMass.h"
#include "ChargeMixCoeff.h"
#include "ChargeMixNdim.h"
#include "ChargeMixRcut.h"
#include "Debug.h"
#include "Dspin.h"
#include "Ecut.h"
#include "Ecutprec.h"
#include "Ecuts.h"
#include "Efield.h"
#include "Polarization.h"
#include "Emass.h"
#include "ExtStress.h"
#include "FermiTemp.h"
#include "IterCmd.h"
#include "IterCmdPeriod.h"
#include "Dt.h"
#include "MuRSH.h"
#include "Nempty.h"
#include "NetCharge.h"
#include "Nrowmax.h"
#include "Nspin.h"
#include "RefCell.h"
#include "ScfTol.h"
#include "Stress.h"
#include "Thermostat.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "Vext.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "Xc.h"

int main(int argc, char **argv, char **envp)
{
  Timer tm;
  tm.start();

#if USE_MPI
  MPI_Init(&argc,&argv);
#endif

  {
  Context ctxt(MPI_COMM_WORLD);

  if ( ctxt.onpe0() )
  {
  cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  cout << "<fpmd:simulation xmlns:fpmd=\"" << qbox_xmlns() << "\">" << endl;
#if USE_UUID
  cout << "<uuid> " << uuid_str() << " </uuid>" << endl;
#endif
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
  cout << "                   I http://qboxcode.org      I\n";
  cout << "                   ============================\n\n";
  cout << "\n";
  cout << "<release> " << release() << " " << TARGET << " </release>" << endl;

  // Identify executable name, checksum, size and link date
  if ( getlogin() != 0 )
    cout << "<user> " << getlogin() << " </user>" << endl;

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
  for ( int i = 0; i < MPI_MAX_PROCESSOR_NAME; i++ )
    processor_name[i] = '\0';
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

#ifdef _OPENMP
  if ( ctxt.onpe0() )
    cout << "<omp_max_threads> " << omp_get_max_threads()
         << " </omp_max_threads>" << endl;
#endif

  UserInterface ui;
  Sample* s = new Sample(ctxt, &ui);

  ui.addCmd(new AngleCmd(s));
  ui.addCmd(new AtomCmd(s));
  ui.addCmd(new BisectionCmd(s));
  ui.addCmd(new ComputeMLWFCmd(s));
  ui.addCmd(new ConstraintCmd(s));
  ui.addCmd(new DistanceCmd(s));
  ui.addCmd(new ExtForceCmd(s));
  ui.addCmd(new FoldInWsCmd(s));
  ui.addCmd(new HelpCmd(s));
  ui.addCmd(new KpointCmd(s));
  ui.addCmd(new ListAtomsCmd(s));
  ui.addCmd(new ListSpeciesCmd(s));
  ui.addCmd(new LoadCmd(s));
  ui.addCmd(new MoveCmd(s));
  ui.addCmd(new PartialChargeCmd(s));
  ui.addCmd(new PlotCmd(s));
  ui.addCmd(new PrintCmd(s));
  ui.addCmd(new QuitCmd(s));
  ui.addCmd(new RandomizeRCmd(s));
  ui.addCmd(new RandomizeVCmd(s));
  ui.addCmd(new RandomizeWfCmd(s));
  ui.addCmd(new RescaleVCmd(s));
  ui.addCmd(new ResetRotationCmd(s));
  ui.addCmd(new ResetVcmCmd(s));
  ui.addCmd(new ResponseCmd(s));
  ui.addCmd(new RseedCmd(s));
  ui.addCmd(new RunCmd(s));
  ui.addCmd(new SaveCmd(s));
  ui.addCmd(new SetCmd(s));
  ui.addCmd(new SetVelocityCmd(s));
  ui.addCmd(new SpeciesCmd(s));
  ui.addCmd(new StatusCmd(s));
  ui.addCmd(new StrainCmd(s));
  ui.addCmd(new TorsionCmd(s));

  ui.addVar(new AlphaPBE0(s));
  ui.addVar(new AlphaRSH(s));
  ui.addVar(new AtomsDyn(s));
  ui.addVar(new BetaRSH(s));
  ui.addVar(new BlHF(s));
  ui.addVar(new BtHF(s));
  ui.addVar(new Cell(s));
  ui.addVar(new CellDyn(s));
  ui.addVar(new CellLock(s));
  ui.addVar(new CellMass(s));
  ui.addVar(new ChargeMixCoeff(s));
  ui.addVar(new ChargeMixNdim(s));
  ui.addVar(new ChargeMixRcut(s));
  ui.addVar(new Debug(s));
  ui.addVar(new Dt(s));
  ui.addVar(new Ecut(s));
  ui.addVar(new Ecutprec(s));
  ui.addVar(new Ecuts(s));
  ui.addVar(new Efield(s));
  ui.addVar(new Polarization(s));
  ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new FermiTemp(s));
  ui.addVar(new IterCmd(s));
  ui.addVar(new IterCmdPeriod(s));
  ui.addVar(new MuRSH(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new NetCharge(s));
  ui.addVar(new Nrowmax(s));
  ui.addVar(new Nspin(s));
  ui.addVar(new Dspin(s));
  ui.addVar(new RefCell(s));
  ui.addVar(new ScfTol(s));
  ui.addVar(new Stress(s));
  ui.addVar(new Thermostat(s));
  ui.addVar(new ThTemp(s));
  ui.addVar(new ThTime(s));
  ui.addVar(new ThWidth(s));
  ui.addVar(new Vext(s));
  ui.addVar(new WfDiag(s));
  ui.addVar(new WfDyn(s));
  ui.addVar(new Xc(s));

  if ( argc == 2 )
  {
    // input file given as a command line argument
    // cmd line: qb inputfilename
    bool echo = true;
    string inputfilename(argv[1]);
    string outputfilename("stdout");
    ifstream in;
    int file_ok = 0;
    if ( ctxt.onpe0() )
    {
      in.open(argv[1],ios::in);
      if ( in )
      {
        // file was opened on process 0
        file_ok = 1;
      }
    }
    MPI_Bcast(&file_ok,1,MPI_INT,0,MPI_COMM_WORLD);

    if ( file_ok )
    {
      ui.processCmds(in, "[qbox]", echo);
    }
    else
    {
      if ( ctxt.onpe0() )
        cout << " Could not open input file " << argv[1] << endl;
    }
  }
  else if ( argc == 4 )
  {
    // server mode
    // cmd line: qb -server inputfilename outputfilename
    if ( strcmp(argv[1],"-server") )
    {
      // first argument is not "-server"
      cout << " use: qb [infile | -server infile outfile]" << endl;
      ctxt.abort(1);
    }
    // first argument is "-server"
    string inputfilename(argv[2]);
    string outputfilename(argv[3]);
    bool echo = false;
    ui.processCmdsServer(inputfilename, outputfilename, "[qbox]", echo);
  }
  else
  {
    // interactive mode
    assert(argc==1);
    // use standard input
    bool echo = !isatty(0);
    string inputfilename("stdin");
    string outputfilename("stdout");
    ui.processCmds(cin, "[qbox]", echo);
  }

  // exit using the quit command when processCmds returns
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( ctxt.onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</fpmd:simulation>" << endl;
  }

  delete s;

  } // end of Context scope
#if USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
