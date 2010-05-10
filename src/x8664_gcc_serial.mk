#-------------------------------------------------------------------------------
#
# Copyright (c) 2009 The Regents of the University of California
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
#  x8664_gcc_serial.mk
#
#-------------------------------------------------------------------------------
# $Id: x8664_gcc_serial.mk,v 1.3 2010-05-10 20:06:42 fgygi Exp $
#
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 MPIDIR=/opt/mpich-1.2.6
 XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_8_0
 FFTWDIR=$(HOME)/software/fftw/Linux_x8664/fftw-2.1.3/fftw
 BLASDIR=/usr/lib64/atlas
 LAPACKDIR=/usr/lib64/atlas

 PLTOBJECTS = readTSC.o

 CXX=/usr/bin/g++
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -g -Wunused -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/.libs -L/usr/X11R6/lib \
           -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lpthread -lfftw \
         -llapack -lf77blas -latlas -lm \
         -Xlinker -Bstatic \
          -lc -lgfortran -static-libgcc -lxerces-c -luuid \
         -Xlinker -Bdynamic

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=Linux_x8664

#-------------------------------------------------------------------------------
