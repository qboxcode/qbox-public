#-------------------------------------------------------------------------------
#
#  aix_mpi.mk
#
#-------------------------------------------------------------------------------
# $Id: aix_mpi.mk,v 1.11 2004-05-20 00:21:18 fgygi Exp $
PLT=AIX
#-------------------------------------------------------------------------------
 XERCESCDIR=${HOME}/software/xml/xerces-c-${PLT}
 XERCESCLIBDIR=/usr/apps/qbox/lib

 CXX=newmpxlC
 LD=$(CXX)

#PLTFLAGS += -DUSE_ESSL -D_LARGE_FILES
 PLTFLAGS += -DUSE_ESSL -DUSE_CSTDIO_LFS
 INCLUDE = -I$(XERCESCDIR)/include
 
 CXXFLAGS= -O2 -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) \
           $(DFLAGS) $(PLTFLAGS)
#CXXFLAGS= -g -qmaxmem=-1 -DUSE_MPI -DSCALAPACK -D$(PLT) $(INCLUDE) $(DFLAGS) \
           $(DFLAGS) $(PLTFLAGS)

 LIBPATH = -L $(XERCESCLIBDIR) 
  
 PLIBS = $(SCALAPACKLIB) $(PBLASLIB) $(TOOLSLIB) $(REDISTLIB) $(CBLACSLIB)
 LIBS =  $(PLIBS) -lessl -lm -lmassv -lxlf90_r $(XERCESCLIBDIR)/libxerces-c.so
 
 LDFLAGS = -bmaxdata:0x80000000 $(LIBPATH) $(LIBS) 

#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = $(HOME)/lib

#
#  MPI setup; uncomment and tailor to your system if using MPIBLACS
#  Will need to comment out the default native BLACS setup below below
#
#USEMPI        = -DUsingMpiBlacs
BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-SP-$(BLACSDBGLVL).a
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-SP-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-SP-$(BLACSDBGLVL).a

#
#  system primitive BLACS setup, comment out if using MPI
#
CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

#
# BLAS and LAPACK
#
#BLASLIB      = /usr/local/lib/lapack.a -lesslsmp
BLASLIB       = /usr/local/lib/lapack.a -lessl

#
#  The name of the libraries to be linked to
#
PBLASLIB      = $(HOME)/lib/pblas_SP.a
SCALAPACKLIB  = $(HOME)/lib/libscalapack_SP.a
TOOLSLIB      = $(HOME)/lib/tools_SP.a
REDISTLIB     = $(HOME)/lib/redist_SP.a
