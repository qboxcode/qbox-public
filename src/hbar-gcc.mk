# ------------------------------------------------------------------------------
#
#  hbar-gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: hbar-gcc.mk,v 1.1 2003-06-04 17:45:39 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/apps/gcc/3.1
 MPIDIR=/home/fgygi/software/gm
 XERCESCDIR=/home/fgygi/software/xml/xerces-c-${PLT}

 CXX=$(GCCDIR)/bin/g++
 LD=$(CXX)

#FFTWDIR=$(HOME)/fftw/linux-pc/fftw-1.3
 FFTWDIR=$(HOME)/fftw/fftw-2.1.3/fftw
 BLASDIR=$(HOME)/software/mkl/lib/32
 
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include
 
 CXXFLAGS= -O3 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) \
            $(INCLUDE) $(DFLAGS) -DAPP_NO_THREADS -DXML_USE_NO_THREADS

# CXXFLAGS= -g  -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS) \
#           -DAPP_NO_THREADS -DXML_USE_NO_THREADS

 LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib -L$(XERCESCDIR)/lib
  
 LIBS =  $(PLIBS) -lfftw -lmkl_lapack $(BLASDIR)/libmkl_def.a \
         -lm -lmpich -lpmpich -lmpich -lgm \
         -lg2c -lxerces-c -lguide -pthread
 
 LDFLAGS = $(LIBPATH) $(LIBS) 

 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = /home/casc/repository/fpmd/software/BLACS/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = /home/casc/repository/fpmd/lib
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/scalapack_$(PLT).a
 TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLT).a
 REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLT).a

 LAPACKLIB = -llapack
 BLASLIB = -lblas

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
