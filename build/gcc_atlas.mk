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
#  gcc_atlas.mk
#
#-------------------------------------------------------------------------------
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 SCALAPACKDIR = $(HOME)/software/scalapack/scalapack-2.0.2
 XERCESCDIR=$(HOME)/software/xerces/xerces-c-src_2_8_0
 ATLASDIR=/usr/lib64/atlas

 PLTOBJECTS = readTSC.o
 SVN_VER :=$(shell svnversion -n)
 DFLAGS += -DSVN_VERSION='"$(SVN_VER)"'

 CXX=mpicxx
 LD=$(CXX)

 PLTFLAGS += -DIA32 -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DMPICH_IGNORE_CXX_SEEK

# FFT must be FFTW2, FFTW3, ESSL or NOLIB
 FFT=FFTW2

ifeq ($(FFT),FFTW2)
 PLTFLAGS += -DUSE_FFTW2
 PLTFLAGS += -DFFTWMEASURE
 FFTWDIR=$(HOME)/software/fftw/Linux_x8664/fftw-2.1.5/fftw
 FFTWINCLUDEDIR=$(FFTWDIR)
 FFTWLIBDIR=$(FFTWDIR)/.libs
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw
endif

ifeq ($(FFT),FFTW3)
 PLTFLAGS += -DUSE_FFTW3
 PLTFLAGS += -DFFTWMEASURE
#PLTFLAGS += -DFFTW_TRANSPOSE
 PLTFLAGS += -DFFTW3_2D
 FFTWDIR=$(HOME)/software/fftw/fftw-3.3.4
 FFTWINCLUDEDIR=$(FFTWDIR)/api
 FFTWLIBDIR=$(FFTWDIR)/.libs
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw3
endif

ifeq ($(FFT),ESSL)
$(error ESSL library not available)
endif

ifeq ($(FFT),NOLIB)
 PLTFLAGS += -DFFT_NOLIB
endif

 INCLUDE += -I$(XERCESCDIR)/include

 CXXFLAGS= -g -O3 -Wunused -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH += -L$(XERCESCDIR)/lib \
            -L$(ATLASDIR) -L$(SCALAPACKDIR) -L/usr/lib64

 LIBS += -lpthread -lxerces-c -lscalapack -llapack -lf77blas -latlas 

 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
