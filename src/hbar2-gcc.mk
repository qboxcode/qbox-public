#------------------------------------------------------------------------------
#
#  hbar2-gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: hbar2-gcc.mk,v 1.1 2005-03-17 23:25:55 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
#GCCDIR=/usr/apps/gcc/gcc-3.4.0
 GMDIR=/GM
 MPIDIR=/MPI
 XERCESCDIR=/home/fgygi/software/xml/xerces-c-src_2_5_0
 FFTWDIR=/usr/apps/fftw
 ATLASDIR=/home/fgygi/software/atlas/Linux_HAMMER64SSE2
 SCALAPACK_DIR = /home/fgygi/software/scalapack/SCALAPACK
 BLACSdir      = /home/fgygi/software/blacs/BLACS/LIB

#CXX=$(GCCDIR)/bin/g++
 CXX=g++
 LD=$(CXX)

 PLTFLAGS = -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
            -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
            -DAPP_NO_THREADS -DXML_USE_NO_THREADS
 
 
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(XERCESCDIR)/include

 CXXFLAGS= -O3 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 

 LIBPATH = -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(ATLASDIR)/lib \
           -L$(XERCESCDIR)/lib \
           -L$(GMDIR)/lib \
#          -L $(GCCDIR)/lib
  
 LIBS =  $(PLIBS) -lfftw -llapack -lf77blas -latlas \
         -lm -lmpich -lpmpich -lmpich -lgm \
         $(XERCESCDIR)/lib/libxerces-c.a  -lg2c
 
 LDFLAGS = $(LIBPATH) $(LIBS) 

 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
# PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/scalapack_$(PLT).a
# TOOLSLIB      = $(SCALAPACK_DIR)/tools_$(PLT).a
# REDISTLIB     = $(SCALAPACK_DIR)/redist_$(PLT).a

#LAPACKLIB = -llapack
#BLASLIB = -lblas
 LAPACKLIB = -llapack
 BLASLIB = -lf77blas

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
