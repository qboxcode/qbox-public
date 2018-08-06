#-------------------------------------------------------------------------------
#
# Copyright (c) 2008-2018 The Regents of the University of California
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
#  centos7_intel.mk
#
#-------------------------------------------------------------------------------
# build on CentOS 7 using the intel compiler and the threaded MKL library
# with openmpi
#
# Prerequisites: 
# On a Centos-7.3 system, install the following packages:
# yum install xerces-c xerces-c-devel
# yum install openmpi openmpi-devel
# yum install fftw fftw-devel
# yum install libuuid libuuid-devel
# intel compiler and MKL library
#
#-------------------------------------------------------------------------------
# Use the following definitions to compile and use qb
# export LD_LIBRARY_PATH=/opt/intel/mkl/intel64/lib
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/lib/intel64
# export PATH=$PATH:/opt/intel/bin

 MKLROOT = /opt/intel/mkl

 PLTOBJECTS = readTSC.o

 CXX=env OMPI_CXX=icc mpicxx
 LD=$(CXX)

 PLTFLAGS += -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES \
             -DXERCESC_3 -DMPICH_IGNORE_CXX_SEEK -DUSE_UUID \
             -DUSE_FFTW3 -DFFTW3_2D

 LIBS = -lfftw3

 INCLUDE = -I$(MKLROOT)/include

 CXXFLAGS= -g -O3 -openmp $(INCLUDE) $(PLTFLAGS) 
 LIBPATH  = -L$(MKLROOT)/lib/intel64

 LIBS +=  -Wl,-Bstatic -Wl,--start-group \
          -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 \
          -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread \
          -Wl,--end-group -Wl,-Bdynamic  \
          -lxerces-c -luuid -liomp5

 LDFLAGS = $(LIBPATH) $(LIBS)

#-------------------------------------------------------------------------------
