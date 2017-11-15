#-------------------------------------------------------------------------------
#
#  bgq.mk
#
#-------------------------------------------------------------------------------
#
 PLT=BGQ
#-------------------------------------------------------------------------------
# use the following .soft environment
####~/.soft###
# @default
# +mpiwrapper-xl 

 CXX=mpic++
 LD=$(CXX)

 PLTFLAGS += -DUSE_MPI -DSCALAPACK
 PLTFLAGS +=  -D__linux__ -DPLT_BIG_ENDIAN
 PLTFLAGS += -DUSE_XERCES 
 PLTFLAGS += -D_LARGEFILE_SOURCE
 PLTFLAGS += -D_FILE_OFFSET_BITS=64
 PLTFLAGS += -DMPICH_IGNORE_CXX_SEEK
 PLTFLAGS += -DUSE_MASSV

# FFT must be FFTW2, FFTW3, ESSL or NOLIB
 FFT=ESSL

ifeq ($(FFT),FFTW2)
 PLTFLAGS += -DUSE_FFTW2
 PLTFLAGS += -DUSE_DFFTW
 PLTFLAGS += -DFFTWMEASURE
 FFTWDIR=/soft/libraries/alcf/current/xl/FFTW2
 FFTWINCLUDEDIR=$(FFTWDIR)/include
 FFTWLIBDIR=$(FFTWDIR)/lib
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -ldfftw
endif

ifeq ($(FFT),FFTW3)
 PLTFLAGS += -DUSE_FFTW3
 PLTFLAGS += -DFFTWMEASURE
 FFTWDIR=/soft/libraries/alcf/current/xl/FFTW3
 FFTWINCLUDEDIR=$(FFTWDIR)/include
 FFTWLIBDIR=$(FFTWDIR)/lib
 INCLUDE += -I$(FFTWINCLUDEDIR)
 LIBPATH += -L$(FFTWLIBDIR)
 LIBS += -lfftw3
endif

ifeq ($(FFT),ESSL)
 PLTFLAGS += -DUSE_ESSL_FFT
 #PLTFLAGS += -DUSE_ESSL_2DFFT
endif

ifeq ($(FFT),NOLIB)
 PLTFLAGS += -DFFT_NOLIB
endif

 XERCESCDIR=$(HOME)/software/xerces/xerces-c-src_2_8_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib

 SCALAPACKLIBDIR=/soft/libraries/alcf/current/xl/SCALAPACK/lib
 SCALAPACKLIB=-lscalapack
 LAPACKLIBDIR=/soft/libraries/alcf/current/xl/LAPACK/lib
 LAPACKLIB=-llapack
 BLASLIBDIR=/soft/libraries/essl/current/lib64
 BLASLIB=-lesslsmpbg

 INCLUDE +=  -I$(XERCESCDIR)/include

 CXXFLAGS= -g -O3 -qsmp=omp -qarch=qp -qtune=qp \
           -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH += -L$(XERCESCLIBDIR) -L$(SCALAPACKLIBDIR) \
            -L$(LAPACKLIBDIR) -L$(BLASLIBDIR) \
            -L${IBM_FCMP}/bglib64 \
            -L${IBM_MAIN_DIR}/xlmass/bg/7.3/bglib64 \
            -L${IBM_MAIN_DIR}/xlf/bg/14.1/bglib64 \
            -L${IBM_MAIN_DIR}/xlsmp/bg/3.1/bglib64 \
            -L$(IBM_MAIN_DIR)/xlf/bg/14.1/bglib64

 LIBS +=  $(PLIBS) $(SCALAPACKLIB) $(LAPACKLIB) $(BLASLIB) \
         -lmass -lmassv  -lxerces-c \
         -lxlf90_r -lxlfmath -lxl \
         -lxlsmp -lxlomp_ser \
         -lpthread -lgomp

 LDFLAGS = $(LIBPATH) $(LIBS)
#-------------------------------------------------------------------------------
