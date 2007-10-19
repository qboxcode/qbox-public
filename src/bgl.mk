#-------------------------------------------------------------------------------
#
#  bgl.mk
#
#-------------------------------------------------------------------------------
# $Id: bgl.mk,v 1.3 2007-10-19 16:24:05 fgygi Exp $
#
 PLT=BGL
#-------------------------------------------------------------------------------
 BGL_ROOT=/bgl/BlueLight/ppcfloor
 BGL_SYS=$(BGL_ROOT)/bglsys

 LIBS_MPI     += -L $(BGL_ROOT)/bglsys/lib -lmpich.rts \
                 -lmsglayer.rts -lrts.rts -ldevices.rts

 GNU_ROOT=/BlueLight/ppcfloor
 BLRTS_GNU_ROOT=$(GNU_ROOT)/blrts-gnu
 CXX=/opt/ibmcmp/vacpp/7.0/bin/blrts_xlC

 LD=$(CXX)
 PLTFLAGS += -DUSE_FFTW \
             -DUSE_MPI -DSCALAPACK \
             -D__linux__ -DPLT_BIG_ENDIAN -DUSE_XERCES \
             -DUSE_CSTDIO_LFS -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64 \
             -DUSE_APC

 FFTWDIR=$(HOME)/software/fftw/bgl/bglfftwgel-2.1.5.pre5
 FFTWINCLUDEDIR=$(FFTWDIR)/fftw
 FFTWLIBDIR=$(FFTWDIR)/fftw/.libs

 APCDIR=$(HOME)/Ctools/counters/bglusr
 APCINCLUDEDIR=$(APCDIR)/include
 APCLIBDIR=$(APCDIR)/lib
 APCLIB=-L$(APCLIBDIR) -lapc.rts -L$(BGL_SYS)/lib -lbgl_perfctr.rts

 XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_6_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib


#BLASDIR=/bgl/local/lib
 BLASDIR=$(HOME)/software/blas/lib

 INCLUDE =  -I$(XERCESCDIR)/include \
            -I$(FFTWINCLUDEDIR) -I$(BGL_ROOT)/bglsys/include \
            -I$(APCINCLUDEDIR)

 CXXFLAGS= -g -O3 -qarch=440 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(FFTWLIBDIR) \
           -L$(BLASDIR) -L$(XERCESCLIBDIR) \
           -L/opt/ibmcmp/xlf/9.1/blrts_lib -L$(APCLIBDIR)

 LIBS =  $(PLIBS) -lfftw $(BLASLIB) -lg2c \
         -lxlf90 -lxlopt -lxlomp_ser -lxl -lxlfmath -lmassv -lxerces-c \
         $(APCLIB)

 LDFLAGS = $(LIBPATH) $(LIBS) $(LIBS_MPI)

 PLAT=BGL
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/bgl/BLACS/LIB
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a
 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/bgl/SCALAPACK
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLAT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

LAPACKLIB = -llapack440
BLASLIB = -ldgemm.rts -lblas440

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
.C.s:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -S $<
#-------------------------------------------------------------------------------
