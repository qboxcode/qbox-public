#-------------------------------------------------------------------------------
#
#  fpmd.mk
#
#-------------------------------------------------------------------------------
# $Id: fpmd.mk,v 1.7 2007-10-19 16:24:05 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/apps/gcc/3.3.2
 MPIDIR=/usr/apps/mpich/1.2.5
 XERCESCDIR=/home/fgygi/software/xml/xerces-c-${PLT}

 CXX=icc
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 PLTOBJECTS = readTSC.o

 FFTWDIR=$(HOME)/fftw/linux-pc-fpmd/fftw-2.1.3/fftw
 BLASDIR=/usr/lib

 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -O3 -xW -Zp16  \
           -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lfftw -llapack -lblas -lm -lmpich -lpmpich -lmpich \
         -lg2c -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS)

 # Blacs libraries
 BLACSDBGLVL   = 0
#BLACSdir      = /home/casc/repository/fpmd/software/BLACS/LIB
 BLACSdir      = $(HOME)/software/blacs/fpmd/BLACS/LIB
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLT)-$(BLACSDBGLVL).a
 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

 # Scalapack libraries
#SCALAPACK_DIR = /home/casc/repository/fpmd/lib
 SCALAPACK_DIR = $(HOME)/software/scalapack/fpmd/SCALAPACK
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

 LAPACKLIB = -llapack
 BLASLIB = -lblas

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
