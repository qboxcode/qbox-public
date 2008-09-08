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
#  mcr-gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: mcr-gcc.mk,v 1.4 2008-09-08 15:56:20 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/local/tools/gnu/gcc/3.2.1_redhat_7a_ia32
 MPIDIR=/usr/lib/mpi
 XERCESCDIR=$(HOME)/software/xml/icc-7.0/xerces-c-src2_2_0

 CXX=$(GCCDIR)/bin/g++
 LD=$(CXX)

 FFTWDIR=$(HOME)/fftw/linux-pc-icc/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/32

 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -O3 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(BLASDIR) -L/usr/lib \
           -L$(XERCESCDIR)/lib

 #LIBS =  $(PLIBS) $(GCCDIR)/lib/libg2c.a -lfftw \
 #        -lmkl_lapack -lmkl -lmkl_def -lmkl_p4 -lm -lmpi -lpmpi \
 #        -lelan -lelan3 -openmp -lrmscall -lxerces-c
 LIBS =  $(PLIBS) $(GCCDIR)/lib/libg2c.a -lfftw \
         -lmkl_lapack -lmkl -lmkl_def -lguide -lpthread -lm -lmpi -lpmpi \
         -lelan -lelan3 -lrmscall -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=INTEL
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/lib
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/lib
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLAT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/scalapack_$(PLAT).a
 TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLAT).a
 REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLAT).a

#LAPACKLIB = -llapack
#BLASLIB = -lblas
 LAPACKLIB = -lmkl_lapack
 BLASLIB = -lmkl

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
