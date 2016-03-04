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
// blacs.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BLACS_H
#define BLACS_H
#include<mpi.h>

extern "C"{
void igesd2d(int*,int*,int*, int*, int*,int*,int*);
void sgesd2d(int*,int*,int*, double*, int*,int*,int*);
void igerv2d(int*,int*,int*, int*, int*,int*,int*);
void sgerv2d(int*,int*,int*, double*, int*,int*,int*);
void sgsum2d(int*,char*,char*,int*,int*,double*,int*,int*,int*);
void igamn2d(int*,char*,char*, int*,
     int*,int*,int*, int*, int*,int*,int*, int*);
void blacs_pinfo(int*, int*);
void blacs_get(int*, int*, int*);
int  sys2blacs_handle(MPI_Comm);
void blacs_barrier(int*, char*);
void blacs_gridinfo(int*, int *, int *, int *, int *);
void blacs_gridinit(int *, char*, int*, int*);
void blacs_gridmap(int*, int *, int*, int*, int*);
void blacs_abort(int*, int*);
void blacs_gridexit(int*);
int  blacs_pnum(int*, int*, int*);
}

#ifdef SCALAPACK
extern "C"{
#endif
// C interface to the BLACS
void Cdgesd2d(int,int,int, double*, int,int,int);
void Cdgerv2d(int,int,int, double*, int,int,int);
void Cdgsum2d(int,char*,char*,int,int,double*,int,int,int);
void Cdgamx2d(int,char*,char*,int,int,double*,int,int*,int*,int,int,int);
void Cdgamn2d(int,char*,char*,int,int,double*,int,int*,int*,int,int,int);
void Cdgebs2d(int,char*,char*,int,int,double*,int);
void Cdgebr2d(int,char*,char*,int,int,double*,int,int,int);

void Cigesd2d(int,int,int, int*, int,int,int);
void Cigerv2d(int,int,int, int*, int,int,int);
void Cigsum2d(int,char*,char*,int,int,int*,int,int,int);
void Cigamx2d(int,char*,char*,int,int,int*,int,int*,int*,int,int,int);
void Cigamn2d(int,char*,char*,int,int,int*,int,int*,int*,int,int,int);
void Cigebs2d(int,char*,char*,int,int,int*,int);
void Cigebr2d(int,char*,char*,int,int,int*,int,int,int);

void Cblacs_pinfo(int*, int*);
void Cblacs_get(int, int, int*);
int  Csys2blacs_handle(MPI_Comm);
void Cblacs_barrier(int, const char*);
void Cblacs_gridinfo(int, int*, int*, int*, int*);
void Cblacs_gridinit(int*, char [], int, int);
void Cblacs_gridmap(int*, int*, int, int, int);
void Cblacs_abort(int, int);
void Cblacs_gridexit(int);
int Cblacs_pnum(int, int, int);
#ifdef SCALAPACK
}
#endif

#endif
