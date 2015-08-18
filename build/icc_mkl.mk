#-------------------------------------------------------------------------------
#
#  icc_mkl.mk
#
#-------------------------------------------------------------------------------
#
 PLT=x86_64
#-------------------------------------------------------------------------------
 MPIDIR=/opt/openmpi
 XERCESCDIR=/share/apps/xerces/xerces-c-src_2_8_0
 PLTOBJECTS = readTSC.o
 INCLUDE = -I$(MKLROOT)/include

 CXX=icc
 LD=mpicxx

 PLTFLAGS += -DIA32 -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

# FFT must be FFTW2, FFTW3, ESSL or NOLIB
  FFT=FFTW2

ifeq ($(FFT),FFTW2)
 PLTFLAGS += -DUSE_FFTW2
 PLTFLAGS += -DFFTWMEASURE
 FFTWDIR=/share/apps/fftw/fftw-2.1.5/fftw
 FFTWINCLUDEDIR=$(FFTWDIR)
 FFTWLIBDIR=$(FFTWDIR)/.libs
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw
endif

ifeq ($(FFT),FFTW3)
 PLTFLAGS += -DUSE_FFTW3
 PLTFLAGS += -DFFTWMEASURE
#PLTFLAGS += -DFFTW_TRANSPOSE
 PLTFLAGS += -DFFTW3_2D
 FFTWDIR=$(HOME)/software/fftw/fftw-3.3.4
 FFTWINCLUDEDIR=$(FFTWDIR)/api
 FFTWLIBDIR=$(FFTWDIR)/.libs
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw3
endif

ifeq ($(FFT),ESSL)
$(error ESSL library not available)
endif

ifeq ($(FFT),NOLIB)
 PLTFLAGS += -DFFT_NOLIB
endif


 INCLUDE += -I$(MPIDIR)/include -I$(XERCESCDIR)/include 

 CXXFLAGS=  -g -O3 -vec-report1 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH += -L$(MPIDIR)/lib64 \
            -L$(MKLDIR)/lib/intel64  \
            -L$(XERCESCDIR)/lib

 LIBS +=  $(PLIBS) \
          -lmkl_intel_lp64 \
          -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core \
          -lirc -lifcore -lsvml \
          -luuid $(XERCESCDIR)/lib/libxerces-c.a -lpthread 

# Parallel libraries
 PLIBS = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
