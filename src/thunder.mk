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
#  thunder.mk
#
#-------------------------------------------------------------------------------
# $Id: thunder.mk,v 1.3 2008-09-08 15:56:20 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/local/tools/gnu/gcc/3.3.3-chaos_2_ia64
 MPIDIR=/usr/lib/mpi
#XERCESCDIR=$(HOME)/software/xml/ia64/xerces-c-src2_2_0-ia64-icc-8.0
 XERCESCDIR=$(HOME)/software/xml/ia64/xerces-c-src_2_5_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib

 CXX=icpc
 LD=$(CXX)

 PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 FFTWDIR=$(HOME)/software/fftw/ia64/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/64

 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -O2 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/.libs \
           -L$(MPIDIR)/lib -L$(BLASDIR) -L/usr/lib \
           -L$(XERCESCLIBDIR)

 LIBS =  $(PLIBS) -lfftw \
         -lmkl -lmkl_lapack -lm -lmpi -lpmpi \
         -lelan -lg2c -openmp $(XERCESCLIBDIR)/libxerces-c.a

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=IA64
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/ia64/BLACS/LIB
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/ia64/SCALAPACK
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

 LAPACKLIB = -lmkl_lapack
 BLASLIB = -lmkl

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
