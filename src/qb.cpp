////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008-2020 The Regents of the University of California
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
// qb.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include <sys/utsname.h>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <cassert>
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
#include "MPIdata.h"

#include "Timer.h"
#include "UserInterface.h"
#include "Sample.h"

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
#include "SpectrumCmd.h"
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
#include "ForceTol.h"
#include "Polarization.h"
#include "Emass.h"
#include "ExtStress.h"
#include "FermiTemp.h"
#include "IterCmd.h"
#include "IterCmdPeriod.h"
#include "LockCm.h"
#include "Dt.h"
#include "MuRSH.h"
#include "Nempty.h"
#include "NetCharge.h"
#include "Nspin.h"
#include "Occ.h"
#include "RefCell.h"
#include "ScfTol.h"
#include "Stress.h"
#include "StressTol.h"
#include "Thermostat.h"
#include "ThTemp.h"
#include "ThTime.h"
#include "ThWidth.h"
#include "Vext.h"
#include "WfDiag.h"
#include "WfDyn.h"
#include "Xc.h"

void usage(void)
{
  cerr << " use:" << endl;
  cerr << "   qb [-nstb nstb] [-nkpb nkpb] [-nspb nspb] inputfile" << endl;
  cerr << "   qb -server [-nstb nstb] [-nkpb nkpb] [-nspb nspb] "
       << "inputfile outputfile" << endl;
}

int main(int argc, char **argv, char **envp)
{
  Timer tm;
  tm.start();

  MPI_Init(&argc,&argv);

  int ntasks;
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  int mype;
  MPI_Comm_rank(MPI_COMM_WORLD,&mype);

  // default values for number of blocks
  // ngb: number of G vector blocks
  // nstb: number of states blocks
  // nkpb: number of kpoint blocks
  // nspb: number of spin blocks
  int ngb = ntasks, nstb = 1, nkpb = 1, nspb = 1;

  // process command line arguments
  int ch;
  opterr = 0; // prevent getopt_long from writing on stderr
  int server_mode = 0;
  bool do_exit = false;

  // options descriptor
  static struct option longopts[] =
  {
    { "help",       no_argument,            NULL,            0  },
    { "server",     no_argument,            NULL,            1  },
    { "nstb",       required_argument,      NULL,            2  },
    { "nkpb",       required_argument,      NULL,            3  },
    { "nspb",       required_argument,      NULL,            4  },
    { NULL,         0,                      NULL,            0  }
  };

  while ( (ch = getopt_long_only(argc, argv, ":", longopts, NULL )) != -1 )
  {
    switch (ch)
    {
      case 0:
        // help
        do_exit = true;
      case 1:
        server_mode = 1;
        break;
      case 2:
        nstb = atoi(optarg);
        if ( nstb < 1 )
        {
          if ( mype == 0 )
            cerr << " nstb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 3:
        nkpb = atoi(optarg);
        if ( nkpb < 1 )
        {
          if ( mype == 0 )
            cerr << " nkpb must be positive" << endl;
          do_exit = true;
        }
        break;
      case 4:
        nspb = atoi(optarg);
        if ( (nspb < 1) || (nspb > 2) )
        {
          if ( mype == 0 )
            cerr << " nspb must be 1 or 2" << endl;
          do_exit = true;
        }
        break;
      case ':':
        if ( mype == 0 )
          cerr << " missing option argument\n";
        do_exit = true;
        break;
      case '?':
        if ( mype == 0 )
          cerr << " unknown or ambiguous option\n";
        do_exit = true;
        break;
      default:
        if ( mype == 0 )
          cerr << " unknown option: " << argv[opterr] << endl;
        do_exit = true;
    }
  }
  argc -= optind;
  argv += optind;

  if ( do_exit )
  {
    if ( mype == 0 )
      usage();
    MPI_Finalize();
    return 0;
  }

  const int interactive = ( argc == 0 );

#ifdef DEBUG
  cout << " argc=" << argc << endl;
  for ( int iarg = 0; iarg < argc; ++iarg )
    cout << " argv[" << iarg << "]= " << argv[iarg] << endl;
  cout << " server_mode = " << server_mode << endl;
  cout << " interactive = " << interactive << endl;
  cout << " ntasks = " << ntasks << endl;
  cout << " nstb = " << nstb << endl;
  cout << " nkpb = " << nkpb << endl;
  cout << " nspb = " << nspb << endl;
#endif

  // adjust ngb to satisfy ngb*nstb*nkpb*nspb == ntasks

  if ( ( ntasks % ( nstb * nkpb * nspb ) ) == 0 )
  {
    // nstb * nkpb * nspb divides ntasks
    ngb = ntasks / ( nstb * nkpb * nspb );
  }
  else
  {
    if ( mype == 0 )
      cerr << " nstb * nkpb * nspb does not divide ntasks evenly" << endl;
    MPI_Finalize();
    return 0;
  }

  MPIdata::set(ngb,nstb,nkpb,nspb);

#ifdef DEBUG
  cout << MPIdata::rank() << ": ngb=" << ngb << " nstb=" << nstb
       << " nkpb=" << nkpb << " nspb=" << nspb << endl;

  cout << MPIdata::rank() << ": igb=" << MPIdata::igb()
       << " istb=" << MPIdata::istb()
       << " ikpb=" << MPIdata::ikpb()
       << " ispb=" << MPIdata::ispb() << endl;
#endif

  if ( MPIdata::onpe0() )
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
    cout << "<release> " << release();
#ifdef TARGET
    cout << " " << TARGET;
#endif
#ifdef VERSION
    cout << " " << VERSION;
#endif
    cout << " </release>" << endl;

    // Identify executable name, checksum, size and link date
    if ( getenv("LOGNAME") != 0 )
      cout << "<user> " << getenv("LOGNAME") << " </user>" << endl;

    // Identify platform
    {
      struct utsname un;
      uname (&un);
      cout << "<sysname> " << un.sysname << " </sysname>" << endl;
      cout << "<nodename> " << un.nodename << " </nodename>" << endl;
    }

    cout << "<start_time> " << isodate() << " </start_time>" << endl;
    cout << " MPIdata::comm: " << ngb << "x" << nstb << "x"
         << nkpb << "x" << nspb << endl;
  }

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

  int coords[4];
  MPI_Cart_coords(MPIdata::comm(),MPIdata::rank(),4,coords);

  if ( MPIdata::onpe0() )
  {
    cout << "<mpi_processes count=\"" << MPIdata::size() << "\">" << endl;
    cout << "<process id=\"" << MPIdata::rank() << "\"> " << processor_name
         << " </process>"
         << " (" << coords[3] << "," << coords[2]
         << "," << coords[1] << "," << coords[0] << ")" << endl;
  }
  for ( int ip = 1; ip < MPIdata::size(); ip++ )
  {
    MPI_Barrier(MPIdata::comm());
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(&buf[0],MPI_MAX_PROCESSOR_NAME,MPI_CHAR,
                   ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(&processor_name[0],MPI_MAX_PROCESSOR_NAME,
        MPI_CHAR,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      MPI_Status status;
      MPI_Recv(coords,4,MPI_INT,ip,ip,MPIdata::comm(),&status);
    }
    else if ( ip == MPIdata::rank() )
    {
      // send processor name to pe0
      MPI_Send(coords,4,MPI_INT,0,MPIdata::rank(),MPIdata::comm());
    }
    if ( MPIdata::onpe0() )
    {
      cout << "<process id=\"" << ip << "\"> " << buf
           << " </process>"
           << " (" << coords[3] << "," << coords[2]
           << "," << coords[1] << "," << coords[0] << ")" << endl;
    }
  }
  if ( MPIdata::onpe0() )
    cout << "</mpi_processes>" << endl;

#ifdef _OPENMP
  if ( MPIdata::onpe0() )
    cout << "<omp_max_threads> " << omp_get_max_threads()
         << " </omp_max_threads>" << endl;
#endif

  UserInterface ui;
  Sample* s = new Sample(&ui);

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
  ui.addCmd(new SpectrumCmd(s));
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
  ui.addVar(new ForceTol(s));
  ui.addVar(new Polarization(s));
  ui.addVar(new Emass(s));
  ui.addVar(new ExtStress(s));
  ui.addVar(new FermiTemp(s));
  ui.addVar(new IterCmd(s));
  ui.addVar(new IterCmdPeriod(s));
  ui.addVar(new LockCm(s));
  ui.addVar(new MuRSH(s));
  ui.addVar(new Nempty(s));
  ui.addVar(new NetCharge(s));
  ui.addVar(new Nspin(s));
  ui.addVar(new Occ(s));
  ui.addVar(new Dspin(s));
  ui.addVar(new RefCell(s));
  ui.addVar(new ScfTol(s));
  ui.addVar(new Stress(s));
  ui.addVar(new StressTol(s));
  ui.addVar(new Thermostat(s));
  ui.addVar(new ThTemp(s));
  ui.addVar(new ThTime(s));
  ui.addVar(new ThWidth(s));
  ui.addVar(new Vext(s));
  ui.addVar(new WfDiag(s));
  ui.addVar(new WfDyn(s));
  ui.addVar(new Xc(s));

  if ( server_mode )
  {
    // server mode
    // input and output files expected as arguments
    if ( argc < 2 )
    {
      if ( MPIdata::onpe0() )
      {
        cout << " server mode requires two arguments" << endl;
        usage();
      }
      MPI_Finalize();
      return 0;
    }
    string inputfilename(argv[0]);
    string outputfilename(argv[1]);
    if ( MPIdata::onpe0() )
    {
      cout << " server mode" << endl;
      cout << " input file:  " << inputfilename << endl;
      cout << " output file: " << outputfilename << endl;
    }
    bool echo = true;
    ui.processCmdsServer(inputfilename, outputfilename, "[qbox]", echo);
  }
  else if ( interactive )
  {
    // interactive mode
    assert(argc==0);
    // use standard input
    bool echo = !isatty(0);
    ui.processCmds(cin, "[qbox]", echo);
  }
  else
  {
    // cmd line: qb inputfilename
    // input file expected as a command line argument
    assert(argc >= 1);
    bool echo = true;
    ifstream in;
    int file_ok = 0;
    if ( MPIdata::onpe0() )
    {
      in.open(argv[0],ios::in);
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
      if ( MPIdata::onpe0() )
        cout << " Could not open input file " << argv[0] << endl;
      MPI_Finalize();
      return 0;
    }
  }

  // exit using the quit command when processCmds returns
  Cmd *c = ui.findCmd("quit");
  c->action(1,NULL);

  if ( MPIdata::onpe0() )
  {
    cout << "<real_time> " << tm.real() << " </real_time>" << endl;
    cout << "<end_time> " << isodate() << " </end_time>" << endl;
    cout << "</fpmd:simulation>" << endl;
  }

  delete s;

  MPI_Finalize();

  return 0;
}
