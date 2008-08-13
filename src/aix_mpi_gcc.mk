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
#  aix_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: aix_mpi_gcc.mk,v 1.3 2008-08-13 06:39:43 fgygi Exp $
PLT=AIX
#-------------------------------------------------------------------------------
 XERCESCDIR=${HOME}/software/xml/xerces-c-${PLT}
 XERCESCLIBDIR=/usr/apps/qbox/lib
 MPIDIR=/usr/local/mpi
 GCCDIR=/usr/local/tools/gnu/gcc/3.1_aix_5

 CXX=$(GCCDIR)/bin/g++
 LD=$(CXX)

 DFLAGS += -DUSE_ESSL -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 INCLUDE = -I$(XERCESCDIR)/include -I$(MPIDIR)/include

#CXXFLAGS= -O2 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)
 CXXFLAGS= -g  -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(GCCDIR)/lib/gcc-lib/powerpc-ibm-aix5.1.0.0/3.1 \
           -L $(XERCESCLIBDIR)  -L $(MPIDIR)/lib

 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)
 LIBS =  $(PLIBS) -lessl -lm -lmassv -lxlf90_r \
         $(XERCESCLIBDIR)/libxerces-c.so -lmpich \
         $(GCCDIR)/lib/libg2c.a

 LDFLAGS = $(LIBPATH) $(LIBS)

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
