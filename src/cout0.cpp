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
// cout0.cpp
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
using namespace std;

#include "MPIdata.h"
#include "cout0.h"

void cout0(string s, int isrc)
{
  // print string s on cout on pe 0 only

  // s contains the msg on rank isrc and is undefined on other ranks
  // isrc is defined on sending rank and is -1 on other ranks
  // send string to pe 0 and print
  bool iamsending = ( isrc == MPIdata::rank() );
  bool iamreceiving = MPIdata::onpe0();
  if ( iamsending && iamreceiving )
  {
    // source == destination: no need to send
    // cout << MPIdata::rank() << ": cout0 len: " << s.size() << endl;
    cout << s;
  }
  else
  {
    // source and destination differ: msg must be sent
    if ( iamreceiving )
    {
      // receive and print
      MPI_Status stat;
      int len;
      MPI_Recv(&len,1,MPI_INT,MPI_ANY_SOURCE,0,MPIdata::comm(),&stat);
      // cout << MPIdata::rank() << ": cout0 len: " << len << endl;
      char* buf = new char[len+1];
      MPI_Recv(buf,len,MPI_CHAR,MPI_ANY_SOURCE,1,MPIdata::comm(),&stat);
      buf[len] = '\0';
      cout << buf;
      delete [] buf;
    }
    if ( iamsending )
    {
      int len = s.size();
      MPI_Send(&len,1,MPI_INT,0,0,MPIdata::comm());
      MPI_Send(s.c_str(),len,MPI_CHAR,0,1,MPIdata::comm());
    }
  }
}
