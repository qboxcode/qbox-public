#-------------------------------------------------------------------------------
#
#  fpmd.mk
#
#-------------------------------------------------------------------------------
# $Id: fpmd.mk,v 1.1 2004-05-20 00:26:18 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/apps/gcc/3.3.2
 MPIDIR=/usr/apps/mpich/1.2.5
 XERCESCDIR=/home/fgygi/software/xml/xerces-c-${PLT}

 CXX=icc
 LD=$(CXX)

 PLTFLAGS += -DUSE_FFTW -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS

 FFTWDIR=$(HOME)/fftw/linux-pc-fpmd/fftw-2.1.3/fftw
 BLASDIR=/usr/lib
 
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include
 
 CXXFLAGS= -O3 -tpp6 -xW -Zp16  \
           -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) 

 LIBPATH = -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L $(BLASDIR) -L $(GCCDIR)/lib -L$(XERCESCDIR)/lib
  
 LIBS =  $(PLIBS) -lfftw -llapack -lblas -lm -lmpich -lpmpich -lmpich \
         -lg2c -lxerces-c
 
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
