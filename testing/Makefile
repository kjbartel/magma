#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2011
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC  = testing_zgemm.cpp	\
	testing_zhemv.cpp	\
	testing_zgemv.cpp	\
	testing_zsymv.cpp	\
	\
	testing_zpotrf.cpp	\
        testing_zpotrf_mgpu.cpp \
	testing_zpotri.cpp	\
	testing_zposv.cpp	\
	testing_zgetrf.cpp 	\
	testing_zgesv.cpp	\
	testing_zgeqrf.cpp 	\
	testing_zgeqlf.cpp 	\
	testing_zgelqf.cpp 	\
	\
	testing_zpotrf_gpu.cpp 	\
	testing_zpotri_gpu.cpp 	\
	testing_zposv_gpu.cpp 	\
	testing_zgetrf_gpu.cpp 	\
        testing_zgetrf_mgpu.cpp \
	testing_zgetri_gpu.cpp 	\
	testing_zgesv_gpu.cpp	\
	testing_zgeqrf_gpu.cpp 	\
        testing_zgeqrf_mgpu.cpp \
        testing_zgelqf_gpu.cpp  \
        testing_zgeqrs_gpu.cpp 	\
        testing_zgeqrs3_gpu.cpp \
	testing_zungqr_gpu.cpp	\
	\
	testing_zgebrd.cpp	\
        testing_zgehrd.cpp 	\
        testing_zhetrd.cpp 	\
        testing_zgeev.cpp       \
        testing_dgeev.cpp       \
        testing_zgesvd.cpp      \
        testing_zheevd.cpp      \
        testing_dsyevd.cpp	\
        testing_zhegvd.cpp      \
        testing_dsygvd.cpp      \
	\
        testing_zhetrd_gpu.cpp  \
        testing_zheevd_gpu.cpp  \
        \
	testing_zcposv_gpu.cpp	\
	testing_zcgesv_gpu.cpp	\
	testing_zcgeqrsv_gpu.cpp

#testing_zswap.cpp

ifeq ($(FORT), pgfortran)
	FOBJ := fortran_thunking.o
	ZSRCF = testing_zgetrf_gpu_f.cuf
else
	FOBJ := fortran.o
	ZSRCF = testing_zgetrf_gpu_f.f90	\
		testing_zgetrf_f.f90
endif

ZSRC += $(ZSRCF)
-include Makefile.src

ALLSRC  = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJF  = $(ALLSRC:.f90=.o)
ALLOBJF := $(ALLOBJF:.cuf=.o)
ALLOBJ   = $(ALLOBJF:.cpp=.o)
EXE      = $(ALLOBJ:.o=)

all: $(EXE)

LIBLAPACKTEST = lin/liblapacktest.a

$(LIBLAPACKTEST):
	cd lin && make liblapacktest.a

$(EXE): $(LIBMAGMA) $(LIBMAGMABLAS) $(LIBLAPACKTEST)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

ifeq ($(FORT), pgfortran)
testing_%.o : testing_%.cuf
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
else
testing_%.o : testing_%.f90
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
endif

%.o: %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

#
# MAGMA Fortran module. Once finalized, the generated magma.mod
# has to be moved in the include directory
#
testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(OPTS) $(LDOPTS) $< $(FOBJ) -o $@ lin/liblapacktest.a -L../lib \
          -lcuda -lmagma -lmagmablas -lmagma $(LIBDIR) $(LIB) 

testing_% : testing_%.o
	$(CC) $(OPTS) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ lin/liblapacktest.a -L../lib \
          -lcuda -lmagma -lmagmablas -lmagma $(LIBDIR) $(LIB) 

#
# These two files are provided by cuda
#
fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
	$(CC) $(OPTS) $(INC) -c $< -o $@
