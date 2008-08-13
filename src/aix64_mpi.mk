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
#  aix64_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: aix64_mpi.mk,v 1.3 2008-08-13 06:39:43 fgygi Exp $
PLT=AIX
#-------------------------------------------------------------------------------
#XERCESCDIR=${HOME}/software/xml/xerces-c-${PLT}
 XERCESCDIR=${HOME}/software/xml/aix64/xerces-c-src2_2_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib

 CXX=mpxlC
 LD=$(CXX)

 PLTFLAGS += -DUSE_ESSL \
             -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 \
             -DUSE_CSTDIO_LFS -DUSE_XERCES -DPLT_BIG_ENDIAN
 INCLUDE = -I$(XERCESCDIR)/include

 CXXFLAGS= -O2 -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) \
           $(DFLAGS) $(PLTFLAGS)
#CXXFLAGS= -g -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS) \
           $(DFLAGS) $(PLTFLAGS)

 LIBPATH = -L $(XERCESCLIBDIR)

#PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)
 LIBS =  $(PLIBS) -lessl -lm -lmassv -lxlf90_r $(XERCESCLIBDIR)/libxerces-c.a

 LDFLAGS = -bmaxdata:0x80000000 $(LIBPATH) $(LIBS)

BLACSDBGLVL   = 0
BLACSdir      = $(HOME)/software/blacs/aix64/BLACS/LIB
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-AIX-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-AIX-$(BLACSDBGLVL).a
CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

BLASLIB       = /usr/local/lib/lapack.a -lessl

SCALAPACKdir  = $(HOME)/software/scalapack/aix64/SCALAPACK
SCALAPACKLIB  = $(SCALAPACKdir)/libscalapack.a
