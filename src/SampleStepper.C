////////////////////////////////////////////////////////////////////////////////
//
// SampleStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SampleStepper.C,v 1.17 2003-11-21 20:01:47 fgygi Exp $

#include "SampleStepper.h"

#include <iostream>
#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SampleStepper::SampleStepper(Sample& s) : s_(s) {}

////////////////////////////////////////////////////////////////////////////////
SampleStepper::~SampleStepper(void)
{
  // print timer map
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( s_.ctxt_.myproc()==0 )
    {
      cout << "<!-- timing "
           << setw(15) << (*i).first
           << " : " << setprecision(3) << setw(9) << tmin
           << " "   << setprecision(3) << setw(9) << tmax << " -->" << endl;
    }
  }
}

