#-------------------------------------------------------------------------------
#
#  mcr.mk
#
#-------------------------------------------------------------------------------
# $Id: mcr_sse2.mk,v 1.1 2004-06-02 21:40:40 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66
 MPIDIR=/usr/lib/mpi
 XERCESCDIR=$(HOME)/software/xml/icc-7.0/xerces-c-src2_2_0
 XERCESCLIBDIR=/usr/apps/qbox/lib
 PLTOBJECTS = readTSC.o

 CXX=icc
 LD=$(CXX)

 PLTFLAGS += -DUSE_SSE2 -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS

 FFTWDIR=$(HOME)/fftw/linux-pc-icc/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/32
 #PAPIDIR=/usr/local/tools/papi
 
 #INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include \
 #          -I$(PAPIDIR)/include
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include 
 
 CXXFLAGS= -g -O3 -xW -Zp16 \
           -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(BLASDIR) -L/usr/lib \
           -L$(XERCESCLIBDIR) 
  
 LIBS =  $(PLIBS) $(GCCDIR)/libg2c.a -lfftw \
         -lmkl_p4 -lmkl_lapack -lm -lmpi -lpmpi \
         -lelan -lelan3 -openmp -lrmscall -lxerces-c
 
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
