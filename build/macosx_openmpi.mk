#
# Copyright (c) 2014 The Regents of the University of California
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
#  macosx_openmpi.mk
#
#-------------------------------------------------------------------------------
#
 PLT=MacOSX_x86_64
#-------------------------------------------------------------------------------
 SCALAPACKDIR=$(HOME)/software/scalapack/scalapack-2.0.2
 LAPACKDIR=$(HOME)/software/lapack/lapack-3.5.0
 BLASDIR=$(HOME)/software/blas/BLAS

 PLTOBJECTS = readTSC.o

 CXX=mpiCC
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_MPI -DUSE_FFTW3 -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES -DXERCESC_3 \
             -DSCALAPACK -DUSE_UUID

 INCLUDE = -I$(FFTWDIR)

 CXXFLAGS= -g -O3 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(SCALAPACKDIR) -L$(LAPACKDIR) -L$(BLASDIR)

 LIBS =  -lfftw3 -lscalapack -llapack -lblas -lm -lgfortran \
         -lxerces-c -lpthread -lstdc++

 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
