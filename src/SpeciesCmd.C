////////////////////////////////////////////////////////////////////////////////
//
// SpeciesCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesCmd.C,v 1.3 2003-05-23 21:51:04 fgygi Exp $

#include "SpeciesCmd.h"
#include "SpeciesReader.h"
using namespace std;

class Species;

////////////////////////////////////////////////////////////////////////////////
int SpeciesCmd::action(int argc, char **argv)
{
  if ( argc != 3 )
  {
    if ( ui->onpe0() )
      cout << "  <!-- use: species name uri -->" << endl;
    return 1;
  }
  
  if ( ui->onpe0() )
    cout << "  <!-- SpeciesCmd: defining species " << argv[1]
         << " as " << argv[2] << " -->" << endl;
         
  SpeciesReader sp_reader(s->ctxt_);
  
  Species* sp = new Species(s->ctxt_,argv[1]);
  
  try
  {
    sp_reader.readSpecies(*sp,argv[2]);
    sp_reader.bcastSpecies(*sp);
    s->atoms.addSpecies(sp,argv[1]);
  }
  catch ( const SpeciesReaderException& e )
  {
    cout << " SpeciesReaderException caught in SpeciesCmd" << endl;
    cout << " SpeciesReaderException: cannot define Species" << endl;
  }
  catch (...)
  {
    cout << " SpeciesCmd: cannot define Species" << endl;
  }
  
  return 0;
}
