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
// SpeciesCmd.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include "SpeciesCmd.h"
#include "SpeciesReader.h"
#include "Species.h"
using namespace std;

class Species;

////////////////////////////////////////////////////////////////////////////////
int SpeciesCmd::action(int argc, char **argv)
{
  if ( argc != 3 )
  {
    if ( ui->onpe0() )
      cout << "  Use: species name uri" << endl;
    return 1;
  }

  if ( ui->onpe0() )
  {
    cout << "  SpeciesCmd: defining species " << argv[1]
         << " as " << argv[2] << endl;
  }

  string xmlstr;
  SpeciesReader sp_reader;
  if ( ui->onpe0() )
    sp_reader.uri_to_string(argv[2], argv[1], xmlstr);

  size_t length = xmlstr.size();
  char* buf = new char[length+1];
  xmlstr.copy(buf,length,0);
  buf[length]='\0';
  MPI_Bcast(buf,xmlstr.size()+1,MPI_CHAR,0,MPI_COMM_WORLD);
  xmlstr = buf;
  delete [] buf;

  //s->ctxt_.string_bcast(xmlstr,0);
  Species* sp = new Species("argv[1]");
  sp_reader.string_to_species(xmlstr,*sp);
  s->atoms.addSpecies(sp,argv[1]);

  return 0;
}
