#-------------------------------------------------------------------------------
#
#  mcr.mk
#
#-------------------------------------------------------------------------------
# $Id: mcr.mk,v 1.12 2007-10-19 16:24:05 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
#GCCDIR=/usr/local/tools/gnu/gcc/3.4.4_chaos_3_x86_elan3/lib
#GCCDIR=/usr/lib/gcc/i386-redhat-linux/3.4.4/libg2c.a
 G2CLIB=-L /usr/local/tools/gnu/gcc/3.4.4_RH_chaos_3_x86/usr/lib -lg2c
 MPIDIR=/usr/lib/mpi
 XERCESCDIR=$(HOME)/software/xml/ia32/chaos_3_x86_elan3/xerces-c-src_2_5_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib
 PLTOBJECTS = readTSC.o

 CXX=icpc
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 FFTWDIR=$(HOME)/software/fftw/ia32/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/32
 #PAPIDIR=/usr/local/tools/papi

 #INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include \
 #          -I$(PAPIDIR)/include
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

#CXXFLAGS= -g -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
 CXXFLAGS= -O3 -xW -Zp16 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(FFTWDIR)/.libs -L/usr/X11R6/lib \
           -L$(BLASDIR) -L$(XERCESCLIBDIR) -L$(MPIDIR)/lib

 LIBS =  $(PLIBS) $(G2CLIB) -lfftw \
        -lmkl_p4 -lmkl_lapack -lm -lmpi -lpmpi \
        -openmp $(XERCESCDIR)/lib/libxerces-c.a

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=IA32
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/ia32/BLACS/LIB
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/ia32/SCALAPACK
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

 LAPACKLIB = -lmkl_lapack
 BLASLIB = -lmkl

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
