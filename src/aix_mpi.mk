#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
#
# This file is distributed under the terms of the GNU General Public
# License. See the file COPYING in the root directory of this distribution
# or http://www.gnu.org/licenses
#
#-------------------------------------------------------------------------------
#
#  aix_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: aix_mpi.mk,v 1.15 2008-06-18 03:39:53 fgygi Exp $
PLT=AIX
#-------------------------------------------------------------------------------
 XERCESCDIR=${HOME}/qb/software/xerces/xerces-c-src_2_7_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib
 FFTWDIR=$(HOME)/qb/software/fftw/fftw-2.1.5/fftw
 FFTWLIBDIR=$(FFTWDIR)/.libs
 FFTWINCLUDEDIR=$(FFTWDIR)

 CXX=mpCC
 LD=$(CXX)

 PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -DUSE_XERCES -DPLT_BIG_ENDIAN
 INCLUDE = -I$(XERCESCDIR)/include -I$(FFTWINCLUDEDIR)

 CXXFLAGS= -O2 -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) \
           $(DFLAGS) $(PLTFLAGS)

 LIBPATH = -L$(XERCESCLIBDIR) -L/usr/local/apps/scalapack -L$(FFTWLIBDIR)

 PLIBS = $(SCALAPACKLIB) $(BLACSLIB)
 LIBS =  $(PLIBS) -lessl -lm -lmassv -lxlf90_r -lfftw \
         $(XERCESCLIBDIR)/libxerces-c.a

 LDFLAGS = -bmaxdata:0x80000000 $(LIBPATH) $(LIBS)

BLACSLIB = -lblacs
SCALAPACKLIB = -lscalapack
