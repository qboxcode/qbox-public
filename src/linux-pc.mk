#-------------------------------------------------------------------------------
#
#  linux-pc.mk
#
#-------------------------------------------------------------------------------
# $Id: linux-pc.mk,v 1.2 2002-10-29 23:46:43 fgygi Exp $
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
