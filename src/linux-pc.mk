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
#  linux-pc.mk
#
#-------------------------------------------------------------------------------
# $Id: linux-pc.mk,v 1.4 2008-08-13 06:39:43 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
 GCCDIR=/usr/apps/gcc/3.1

 CXX=$(GCCDIR)/bin/g++
 LD=$(CXX)

 FFTWDIR=$(HOME)/fftw/linux-pc/fftw-1.3/src
 BLASDIR=/usr/lib

 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR)

 CXXFLAGS= -O2 -D$(PLT) -DADD_ $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(BLASDIR) -L$(FFTWDIR) -L$(GCCDIR)/lib
#LIBPATH = -L$(BLASDIR) -L$(FFTWDIR) -L/usr/X11R6/lib -L$(GCCDIR)/lib

 LIBS = -lfftw -llapack -lblas -lm -lg2c

 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
