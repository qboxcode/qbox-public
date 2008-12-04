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
#  x8664_gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: x8664_gcc.mk,v 1.14 2008-12-04 20:05:13 fgygi Exp $
#
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 MPIDIR=/opt/mpich-1.2.6
#XERCESCDIR=$(HOME)/software/xml/Linux_x8664/xerces-c-src_2_5_0
 XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_8_0
 FFTWDIR=$(HOME)/software/fftw/Linux_x8664/fftw-2.1.3/fftw
 BLASDIR=$(HOME)/software/atlas/ATLAS/Linux_P4E64SSE3/lib
 LAPACKDIR=$(HOME)/software/lapack/LAPACK

 PLTOBJECTS = readTSC.o

 CXX=/usr/bin/g++
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -g -Wunused -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/.libs -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lpthread -lfftw \
         -llapack -lf77blas -latlas -lm \
         -Xlinker -Bstatic \
          -lc -lgfortran -static-libgcc -lmpich -lxerces-c \
         -Xlinker -Bdynamic

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/Linux_x8664/BLACS/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/Linux_x8664/scalapack-1.8.0
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
