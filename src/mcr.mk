#-------------------------------------------------------------------------------
#
#  linux-pc_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: mcr.mk,v 1.2 2003-11-21 19:09:01 fgygi Exp $
#
 PLT=LINUX
#-------------------------------------------------------------------------------
 GCCDIR=/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66
 MPIDIR=/usr/lib/mpi
 XERCESCDIR=$(HOME)/software/xml/icc-7.0/xerces-c-src2_2_0
 XERCESCLIBDIR=/usr/apps/qbox/lib

 CXX=icc
 LD=$(CXX)

 DFLAGS += -DUSE_CSTDIO_LFS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
 
 FFTWDIR=$(HOME)/fftw/linux-pc-icc/fftw-2.1.3/fftw
 BLASDIR=/opt/intel/mkl/lib/32
 
 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include
 
#CXXFLAGS= -O2 -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)
 CXXFLAGS= -O3 -xW -Zp16 -tpp7  -DUSE_MPI -DSCALAPACK -DADD_ -D$(PLT) $(INCLUDE) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(BLASDIR) -L/usr/lib \
           -L$(XERCESCLIBDIR)
  
 LIBS =  $(PLIBS) $(GCCDIR)/libg2c.a -lfftw \
         -lmkl_lapack -lmkl -lmkl_def -lmkl_p4 -lm -lmpi -lpmpi \
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
