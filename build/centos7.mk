#-------------------------------------------------------------------------------
#
# Copyright (c) 2008-2017 The Regents of the University of California
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
#  centos7.mk
#
#-------------------------------------------------------------------------------
# Prerequisites: 
# On a Centos-7.3 system, install the following packages:
# yum install xerces-c xerces-c-devel
# yum install openmpi openmpi-devel
# yum install lapack lapack-devel
# yum install fftw fftw-devel
# yum install scalapack-common scalapack-openmpi \
#             scalapack-openmpi-devel scalapack-openmpi-static
# yum install libuuid libuuid-devel
#
#-------------------------------------------------------------------------------
 PLT=Linux_x8664
#-------------------------------------------------------------------------------

 PLTOBJECTS = readTSC.o

 CXX=mpicxx
 LD=$(CXX)

 PLTFLAGS += -DIA32 -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DXERCESC_3 -DMPICH_IGNORE_CXX_SEEK -DUSE_UUID

# FFT must be FFTW2, FFTW3, ESSL or NOLIB
 FFT=FFTW3

ifeq ($(FFT),FFTW2)
 PLTFLAGS += -DUSE_FFTW2
#PLTFLAGS += -DFFTWMEASURE
#FFTWDIR=$(HOME)/software/fftw/Linux_x8664/fftw-2.1.5/fftw
#FFTWINCLUDEDIR=$(FFTWDIR)
#FFTWLIBDIR=$(FFTWDIR)/.libs
#INCLUDE += -I$(FFTWINCLUDEDIR)
#LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw
endif

ifeq ($(FFT),FFTW3)
 PLTFLAGS += -DUSE_FFTW3
#PLTFLAGS += -DFFTWMEASURE
#PLTFLAGS += -DFFTW_TRANSPOSE
 PLTFLAGS += -DFFTW3_2D
#FFTWDIR=$(HOME)/software/fftw/fftw-3.3.4
#FFTWINCLUDEDIR=$(FFTWDIR)/api
#FFTWLIBDIR=$(FFTWDIR)/.libs
#INCLUDE += -I$(FFTWINCLUDEDIR)
#LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw3
endif

ifeq ($(FFT),ESSL)
$(error ESSL library not available)
endif

ifeq ($(FFT),NOLIB)
 PLTFLAGS += -DFFT_NOLIB
endif

 CXXFLAGS= -g -O3 -Wunused -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
 LIBS += -lpthread -lxerces-c -lscalapack -llapack -lblas -luuid
 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
