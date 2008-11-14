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
# $Id: aix_mpi.mk,v 1.18 2008-11-14 04:03:20 fgygi Exp $
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
