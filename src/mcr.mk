#-------------------------------------------------------------------------------
#
#  mcr.mk
#
#-------------------------------------------------------------------------------
# $Id: mcr.mk,v 1.9 2004-09-14 22:24:11 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/lib/gcc-lib/i386-redhat-linux/3.2.3
 MPIDIR=/usr/lib/mpi
 XERCESCDIR=$(HOME)/software/xml/xerces-c-src2_2_0
#XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_5_0
 XERCESCLIBDIR=/usr/apps/qbox/lib
 PLTOBJECTS = readTSC.o

 CXX=icc
 LD=$(CXX)

 PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 FFTWDIR=$(HOME)/fftw/linux-pc-icc/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/32
 #PAPIDIR=/usr/local/tools/papi
 
 #INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include \
 #          -I$(PAPIDIR)/include
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include 
 
#CXXFLAGS= -g -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 
 CXXFLAGS= -O3 -xW -Zp16 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(BLASDIR) \
           -L$(XERCESCLIBDIR) 
  
 LIBS =  $(PLIBS) $(GCCDIR)/libg2c.a -lfftw \
        -lmkl_p4 -lmkl_lapack -lm -lmpi -lpmpi \
        -openmp $(XERCESCDIR)/lib/libxerces-c.a
 
 LDFLAGS = $(LIBPATH) $(LIBS) 

 PLAT=INTEL
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/lib
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/lib
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLAT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/scalapack_$(PLAT).a
 TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLAT).a
 REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLAT).a

#LAPACKLIB = -llapack
#BLASLIB = -lblas
 LAPACKLIB = -lmkl_lapack
 BLASLIB = -lmkl

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
