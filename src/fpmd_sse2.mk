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
#  fpmd.mk
#
#-------------------------------------------------------------------------------
# $Id: fpmd_sse2.mk,v 1.4 2008-09-08 15:56:19 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/apps/gcc/3.3.2
 MPIDIR=/usr/apps/mpich/1.2.5
 XERCESCDIR=/home/fgygi/software/xml/xerces-c-${PLT}

 CXX=icc
 LD=$(CXX)

 PLTFLAGS += -DUSE_SSE2 -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS

 PLTOBJECTS = readTSC.o

 FFTWDIR=$(HOME)/fftw/linux-pc-fpmd/fftw-2.1.3/fftw
 BLASDIR=/usr/lib

 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -O3 -xW -Zp16  \
           -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lfftw -llapack -lblas -lm -lmpich -lpmpich -lmpich \
         -lg2c -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS)

 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = /home/casc/repository/fpmd/software/BLACS/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = /home/casc/repository/fpmd/lib
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/scalapack_$(PLT).a
 TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLT).a
 REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLT).a

 LAPACKLIB = -llapack
 BLASLIB = -lblas

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
