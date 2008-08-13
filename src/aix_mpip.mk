#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
#
# This file is part of Qbox
#
# Qbox is distributed under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
# See the file COPYING in the root directory of this distribution
# or <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
#
#  aix_mpip.mk
#
#-------------------------------------------------------------------------------
# $Id: aix_mpip.mk,v 1.5 2008-08-13 06:39:43 fgygi Exp $
PLT=AIX
#-------------------------------------------------------------------------------
 XERCESCDIR=${HOME}/software/xml/xerces-c-${PLT}

 CXX=newmpxlC
 LD=$(CXX)

 INCLUDE = -I$(XERCESCDIR)/include

 CXXFLAGS= -O3 -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
#CXXFLAGS= -g -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(XERCESCDIR)/lib  -L/usr/local/mpiP/lib

 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)
 LIBS =  $(PLIBS) -lessl -lm -lmassv \
         -lxlf90_r -lmpiP -lbfd -liberty -lintl \
         ${XERCESCDIR}/lib/libxerces-c.so

 LDFLAGS = -bmaxdata:0x80000000 $(LIBPATH) $(LIBS)

#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = $(HOME)/lib

#
#  MPI setup; uncomment and tailor to your system if using MPIBLACS
#  Will need to comment out the default native BLACS setup below below
#
#USEMPI        = -DUsingMpiBlacs
BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-SP-$(BLACSDBGLVL).a
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-SP-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-SP-$(BLACSDBGLVL).a

#
#  system primitive BLACS setup, comment out if using MPI
#
CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

#
# BLAS and LAPACK
#
#BLASLIB      = /usr/local/lib/lapack.a -lesslsmp
BLASLIB       = /usr/local/lib/lapack.a -lessl

#
#  The name of the libraries to be linked to
#
PBLASLIB      = $(HOME)/lib/pblas_SP.a
SCALAPACKLIB  = $(HOME)/lib/libscalapack_SP.a
TOOLSLIB      = $(HOME)/lib/tools_SP.a
REDISTLIB     = $(HOME)/lib/redist_SP.a
